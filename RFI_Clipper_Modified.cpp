/**
*********************************************************************
1. Include predefined libraries
*********************************************************************
**/
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iterator>
#include <algorithm>
#include "omp.h"
#include <iostream>
#include <vector>
#define BOOST_CB_DISABLE_DEBUG // The Debug Support has to be disabled, otherwise the code produces a runtime error.
#include <boost/circular_buffer.hpp>
#include <boost/assert.hpp>
#include <assert.h>
#include <ctime>
using namespace std;


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}


/**
*********************************************************************
2. Beginning of the main function.
*********************************************************************
**/
int main(int argc, char * argv [])
{
  
  char * infilename = getCmdOption(argv, argv + argc, "-inputFile");
  char * outfilename = getCmdOption(argv, argv + argc, "-outputFile");
  char * directory = getCmdOption(argv, argv + argc, "-outputDirectory");
  char * fileNumber = getCmdOption(argv, argv + argc, "-fileNumber");


   time_t tstart, tend; 
   tstart = time(0);

  /**
  *********************************************************************
  3. Set the default values before processing
  *********************************************************************
  **/
  int nSubbands   = 1;
  unsigned goodSamples = 0;
  float bandPassAve = 110.0;
  float bandPassVar = 16.0;
  float dataChunksInSec = 10;
  float spectrumRMS = 0.0;
  float spectrumSum = 0.0;
  float dataModel   = 0.0;
  float _rmsRunAve  = 0.0;
  float _meanRunAve = 0.0;
  float _useMeanOverRMS = 0;
  float  _crFactor = 10; // channel rejection RMS 
  float _srFactor = 4; // spectrum rejection RMS
  float _zeroDMing = 1;
  float k = 4; // degrees of freedom
  int nbands = 8;
  float meanMinusMinimum = k / sqrt(2.0*k);
  std::vector<float> _lastGoodSpectrum;
  int _badSpectra(0);
  int _current(0); // history pointer
  int NumberOfFlaggedSamples = 0;
  std::vector<float> SamplesReplacedType1, SamplesReplacedType2, SamplesReplacedType3;
  boost::circular_buffer<float> _meanBuffer, _rmsBuffer;
  _meanBuffer.set_capacity(_maxHistory);
  _rmsBuffer.set_capacity(_maxHistory);
  std::vector<unsigned char> n(1000);
  std::vector<unsigned char> tsampRaw(8);
  std::vector<unsigned char> nbitsRaw(4);
  std::vector<unsigned char> nchansRaw(4);
  double tsamp,_meanOverRMS, _lastGoodMean, _lastGoodRMS;
  float _fractionBadChannels;
  int nChannels, nbits, headerSize, plek, _remainingZeros;
  bool _num = 0;
  /**
  ****************************************************************************
  4. Read the header of the SIGPROC file and write the header to the new file
  ****************************************************************************
  **/

  // Determine the size of the whole file
  std::streampos fsize = 0;
  ifstream infile;
  infile.open(infilename, ios::binary | ios::in);
  //infile.open("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0.fil", ios::binary | ios::in);
  fsize = infile.tellg();
  infile.seekg( 0, std::ios::end );
  fsize = infile.tellg() - fsize;
  std::cout << "The size of the file is: " << fsize << std::endl;

  //Determine the size of the header
  infile.seekg(0, ios::beg);
  infile.read(reinterpret_cast<char*>(n.data()), n.size()*sizeof(unsigned char));
  const string headerEND ="END";
  auto idx1 = std::search(n.begin(), n.end(), headerEND.begin(), headerEND.end());
  headerSize = idx1-n.begin()+3;
  std::cout<< "Size of the header in bytes: "<< headerSize << std::endl;
  
  //Determine the number of bits
  const string StringNbits ="nbits";
  auto idx2 = std::search(n.begin(), n.end(), StringNbits.begin(), StringNbits.end());
  plek = idx2-n.begin()+5;
  infile.seekg(plek, ios::beg);
  infile.read(reinterpret_cast<char*>(nbitsRaw.data()), nbitsRaw.size()*sizeof(unsigned char));
  std::copy(reinterpret_cast<const char*>(&nbitsRaw[0]), reinterpret_cast<const char*>(&nbitsRaw[4]), reinterpret_cast<char*>(&nbits));
  std::cout<< "Number of bits are: "<< nbits << std::endl;

  //Write the header to the new file 
  ofstream outfile;
  outfile.open(outfilename, ios::binary);
  //outfile.open("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_normalised.fil", ios::binary);
  std::vector<unsigned char> header(headerSize);
  infile.seekg(0, ios::beg);
  infile.read(reinterpret_cast<char*>(header.data()), header.size()*sizeof(unsigned char));
  outfile.write(reinterpret_cast<char*>(header.data()), header.size()*sizeof(unsigned char));
  int normalisedDataSavedAsFloat = 32;
  outfile.seekp(plek, ios::beg);
  outfile.write(reinterpret_cast<const char *>(&normalisedDataSavedAsFloat), sizeof(normalisedDataSavedAsFloat));
  outfile.close();

  //Determine the number of channels
  const string StringNchans="nchans";
  auto idx3 = std::search(n.begin(), n.end(), StringNchans.begin(), StringNchans.end());
  plek = idx3-n.begin()+6;
  infile.seekg(plek, ios::beg);
  infile.read(reinterpret_cast<char*>(nchansRaw.data()), nchansRaw.size()*sizeof(unsigned char));
  std::copy(reinterpret_cast<const char*>(&nchansRaw[0]), reinterpret_cast<const char*>(&nchansRaw[4]), reinterpret_cast<char*>(&nChannels));
  unsigned nBins = nChannels * nSubbands;
  std::cout<< "Number of channels are: "<< nChannels <<"\n";

  //Determine the sampling interval
  const string TsampInterval ="tsamp";
  auto idx4 = std::search(n.begin(), n.end(), TsampInterval.begin(), TsampInterval.end());
  plek = idx4-n.begin()+5;
  infile.seekg(plek, ios::beg);
  infile.read(reinterpret_cast<char*>(tsampRaw.data()), tsampRaw.size()*sizeof(unsigned char));
  std::copy(reinterpret_cast<const char*>(&tsampRaw[0]), reinterpret_cast<const char*>(&tsampRaw[8]), reinterpret_cast<char*>(&tsamp));
  std::cout<< "The sampling interval is: "<< tsamp<<"\n";
  infile.close();
  //nSamples = (float(fsize)-headerSize)/(nbits/8.0)/nChannels-1;
  //std::cout<<"Number of samples are: "<<nSamples<<std::endl;
  //nSamples = 100;

  float _maxHistory  = 1/tsamp;  // Average filter length
  int tobs = (float(fsize)-headerSize)/(nChannels*nSubbands)*tsamp;
  int nSamples = (float(fsize)-headerSize)/(nbits/8.0)/nChannels-1;
  int NumberOfChunks  = floor(tobs/dataChunksInSec);
  int NumberOfSamplesPerChunk = floor(dataChunksInSec/tsamp);
  int averageSamplestoReadFromFile = (nChannels*NumberOfSamplesPerChunk);


  std::string nameofoutput2 =  directory + std::string(fileNumber) + std::string("_RFI_Type2.txt");
  std::string nameofoutput1 =  directory + std::string(fileNumber) + std::string("_RFI_Type1.txt"); 
  std::string nameofoutput3 =  directory + std::string(fileNumber) + std::string("_RFI_Type3.txt"); 

  for (unsigned a=0 ; a < NumberOfChunks ; ++a)
  {  
    /**
    ****************************************************************************
    5. Read in the data to be processed. Data is stored in "I"
    ****************************************************************************
    **/
    infile.open(infilename, ios::binary | ios::in);
    std::vector<unsigned char> dataRaw(averageSamplestoReadFromFile);
    infile.seekg(headerSize+a*averageSamplestoReadFromFile, ios::beg);
    infile.read(reinterpret_cast<char*>(dataRaw.data()), dataRaw.size()*sizeof(unsigned char));
    std::vector<float> I(dataRaw.begin(), dataRaw.end());
    dataRaw.clear();
    infile.close();

    if (_lastGoodSpectrum.size() != nBins)
    {
      _lastGoodSpectrum.resize(nBins,0.0);
      _remainingZeros = nBins;
      std::cout << "RFI_Clipper: resizing _lastGoodSpectrum" << std::endl;
    }

    /**
    ****************************************************************************
    6. Start processing the data: main FOR-loop starts here.
    ****************************************************************************
    **/
    for (unsigned t = 0; t < NumberOfSamplesPerChunk; ++t) 
    {
      const std::vector<float> bandPass(nChannels,bandPassAve);
      float spectrumSum = 0.0;
      float spectrumSumSq = 0.0;
      // Split the spectrum into 8 bands for the purpose of matching
      // the spectrum to the model
      std::vector<float> goodChannels(8,0.0);
      unsigned channelsPerBand = nBins / 8;
      // find the minima of I in each band, and the minima of the
      // model and compare
      std::vector<float> miniData(8,1e6);
      std::vector<float> miniModel(8,1e6);
      std::vector<float> dataMinusModel(8);
      std::vector<float> bandSigma(8);
      std::vector<float> bandMean(8,0.0);
      std::vector<float> bandMeanSquare(8,0.0);


      // Find the data minima and model minima in each band
      // Let us also estimate sigma in each band

      // Try this over an adapting stage, lasting as long as the running average buffers
      if (!_rmsBuffer.full())
      {
        for (unsigned s = 0; s < nSubbands; ++s) 
        {
          long index = t*nChannels*nSubbands+s*nChannels;
          for (unsigned c = 0; c < nChannels ; ++c) 
          {
            int binLocal = s*nChannels +c;
            unsigned band = (int)binLocal / channelsPerBand;

            //    std::cout << "----------------" << std::endl; std::
            //    cout << band << " " << binLocal << " " << I[index+c]
            //    << std::endl;

            if (I[index+c] < miniData[band]) miniData[band] = I[index+c];
            if (bandPass[binLocal] < miniModel[band]) miniModel[band] = bandPass[binLocal];
            bandMean[band] += I[index+c];
            bandMeanSquare[band] += (I[index+c]*I[index+c]);
          }
        }


        // Now find the distances between data and model and the RMSs in
        // each band

        for (unsigned b = 0; b < 8; ++b)
        {
          dataMinusModel[b] = miniData[b] - miniModel[b];
          bandSigma[b] = sqrt(bandMeanSquare[b]/channelsPerBand - std::pow(bandMean[b]/channelsPerBand,2));
          //std::cout << bandSigma[b] << " " << bandMean[b] << " " <<
          //bandMeanSquare[b] << std::endl;
        }

        // Assume the minimum bandSigma to be the best estimate of this
        // spectrum RMS
        spectrumRMS = *std::min_element(bandSigma.begin(), bandSigma.end());
        
        // Take the median of dataMinusModel to determine the distance from the model
        std::nth_element(dataMinusModel.begin(), dataMinusModel.begin()+dataMinusModel.size()/2, dataMinusModel.end());
        dataModel = (float)*(dataMinusModel.begin()+dataMinusModel.size()/2);
        // std::cout << "data minus model " << dataModel << " " <<
        // std::endl; since we have used the minima to determine this
        // distance, we assume that dataModel is actually k/sqrt(k)
        // sigma away from the real value, where k is the number of the
        // degrees of freedom of the chi-squared distribution of the
        // incoming data. For no integration, k will be 4 (2 powers per poln)

        // Let us now build up a running average of spectrumRMS values
        // (_maxHistory of them)

        // if the buffer is not full, compute the new rmsRunAve like this
        // if (!_rmsBuffer.full()) {
        _rmsBuffer.push_back(spectrumRMS);
        _rmsRunAve = std::accumulate(_rmsBuffer.begin(), _rmsBuffer.end(), 0.0)/_rmsBuffer.size();

        // The very last time this is done, store a reference value of
        // the mean over the rms; this works as I have just added on the
        // last value two lines above
        if (_rmsBuffer.full()) 
        {
          _meanOverRMS = _meanRunAve / _rmsRunAve;
          dataModel = 0;
        }

        // and update the model
        dataModel = dataModel + meanMinusMinimum * _rmsRunAve;
      }
      else
      {
        // just update the running average with the current last value,
        // and take the oldest off the end. Then add the new
        // value onto the buffer. The idea here is that the new
        // value is not used for the current spectrum, but rather
        // for the one after. The current spectrum is evaluated
        // based on the rms values in the buffer up to that
        // point.

        _rmsRunAve -= (_rmsBuffer.front()/_rmsBuffer.size());
        _rmsRunAve += (_rmsBuffer.back()/_rmsBuffer.capacity());

        // In extreme RFI cases, the measured RMS may start growing
        // due to particular RFI signals. The mean over rms ratio
        // should remain approximately constant over the course of the
        // observation. Use this as a safety check, after the mean has
        // been reasonably determined, to set the RMS to a more
        // reliable value :
        if  (_useMeanOverRMS)
           //recover the rms running average
          spectrumRMS = std::abs(_meanRunAve / _meanOverRMS);

        _rmsBuffer.push_back(spectrumRMS);
        // use the mean running average as a model of the mean of the
        // data; remember that the running average is updated with a
        // mean after channel clipping, below.
        //dataModel = _meanRunAve;
      }

      // now use this rms to define a margin of tolerance for bright
      // channels
      float margin = _crFactor * _rmsRunAve;

      // Now loop around all the channels: if you find a channel where
      // (I - bandpass) - datamodel > margin, then replace it
      for (unsigned s = 0; s < nSubbands; ++s) 
      {
        long index = t*nChannels*nSubbands+s*nChannels;
        for (unsigned c = 0; c < nChannels; ++c) 
        {
          int binLocal = s*nChannels +c;
          //std::cout<< "Raw data: "<< index+c <<" "<<I[index+c] <<std::endl;
          if (I[index+c] - dataModel - bandPass[binLocal] > margin) 
          {
            //std::ofstream f1("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type1.txt", ios::app);
            std::ofstream f1(nameofoutput1, ios::app);
            // clipping this channel to values from the last good
            // spectrum
            /*std::cout << "clipping channel " << I[index+c]
              << " " << dataModel << " " <<  bandPass[binLocal]
              << " " << margin << std::endl;*/
            //The following is for polarization
            I[index + c] = _lastGoodSpectrum[binLocal];
            long place = a*averageSamplestoReadFromFile+index + c ;
            f1 << (place) << "\n";
            //SamplesReplacedType1.push_back(place);
            f1.close();
          }
          else
          {
            unsigned band = (int)binLocal / channelsPerBand;
            ++goodChannels[band];
            spectrumSum += I[index+c];
          }
        }
      }
      // So now we have the mean of the incoming data, in a reliable
      // form after channel clipping
      unsigned totalGoodChannels=std::accumulate(goodChannels.begin(), goodChannels.end(), 0);
      _fractionBadChannels += (float)(nBins - totalGoodChannels)/nBins;
      spectrumSum /= totalGoodChannels;

      // Check if more than 20% of the channels in each band were
      // bad. If so in more than half of the bands, 4 in this case,
      // keep record. Also, if one band is completely gone, or less
      // than 80% of the total survive, keep record.
      unsigned badBands = 0;
      for (unsigned b = 0; b < 8; ++b)
      {
        if (goodChannels[b] < 0.8 * channelsPerBand) 
        {
          ++badBands;
        }
        if (goodChannels[b] == 0) 
        {
          badBands += 4;
        }
      }
      if (totalGoodChannels < 0.8 * nBins) badBands +=4;

      // Let us now build up the running average of spectrumSum values
      // (_maxHistory of them) if the buffer is not full, compute the
      // new meanRunAve like this
      if (!_meanBuffer.full()) 
      {
        _meanBuffer.push_back(spectrumSum);
        _meanRunAve = std::accumulate(_meanBuffer.begin(), _meanBuffer.end(), 0.0)/_meanBuffer.size();
      }
      else 
      {
        //   just update the running average with the new value, and
        //   take the oldest off the end, using the same principle as
        //   with the rms buffer, i.e. do not use the current
        //   measurement for the current spectrum.

        // Note there is a tiny descrepance at the point when the
        // buffer is first full

        //  std::cout << "History buffer now full " << _num << std::endl;
        _meanRunAve -= _meanBuffer.front()/_meanBuffer.size();
        _meanRunAve += _meanBuffer.back()/_meanBuffer.size();
        _meanBuffer.push_back(spectrumSum);
      }

      // Now we can check if this spectrum has an abnormally high mean
      // compared to the running average

      // Let us define the tolerance first, and remember, we are
      // summing across nBins, hence sqrt(nBins), strictly only valid
      // for Gaussian stats
      float spectrumRMStolerance = _srFactor * bandPassVar/sqrt(nBins);

      //Now check, if spectrumSum - model > tolerance, declare this
      //time sample useless, replace its data and take care of the
      //running averages, also cut spectra where badBands >= 4, see above

      if (spectrumSum - _meanRunAve > spectrumRMStolerance || badBands >= 4) 
      {
        //std::string nameofoutput2 =  "/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SN" + std::to_string(0) + "_RFI_Type2.txt";
        std::ofstream f2(nameofoutput2, ios::app);
        // we need to remove this entire spectrum
        for (unsigned s = 0; s < nSubbands; ++s) 
          {
            long index = t*nChannels*nSubbands+s*nChannels;
            for (unsigned c = 0; c < nChannels; ++c) 
            {
              int binLocal = s*nChannels +c;
              long place = a*averageSamplestoReadFromFile+index + c ;
              I[index + c] = _lastGoodSpectrum[binLocal];
              f2 << (place) << "\n";
              //cout << place << std::endl;
              //SamplesReplacedType2.push_back((place));
            }
          }
        f2.close();
        // keep a record of bad spectra
        ++_badSpectra;

        // now remove the last samples from the running average
        // buffers and replace them with the last good values.
        //
        _meanBuffer.pop_back();
        _rmsBuffer.pop_back();
        //  _meanBuffer.push_back(_lastGoodMean);
        //  _rmsBuffer.push_back(_lastGoodRMS);
        _meanBuffer.push_back(_meanBuffer.back());
        _rmsBuffer.push_back(_rmsBuffer.back());
      }
      // else keep a copy of the original spectrum, as it is good
      else
      {
        spectrumSum = 0.0;
        spectrumSumSq = 0.0;
        //  if (!_lastGoodSpectrum.full())
        // if _lastGoodSpectrum is full, _remainingZeros will be 0
        if (_remainingZeros != 0)
        {
          for (unsigned s = 0; s < nSubbands; ++s) 
          {
            long index = t*nChannels*nSubbands+s*nChannels;
            for (unsigned c = 0; c < nChannels; ++c) 
            {
              int binLocal = s*nChannels +c;
              if (_lastGoodSpectrum[binLocal] == 0.0 && I[index+c] != 0.0 )
              {
                _lastGoodSpectrum[binLocal] = I[index+c];
                --_remainingZeros;
                std::cout << "remaining channels to fill in good spectrum: "
                << _remainingZeros << " " << _num
                << " " << totalGoodChannels << std::endl;
                spectrumSum += _lastGoodSpectrum[binLocal];
                spectrumSumSq += _lastGoodSpectrum[binLocal] *_lastGoodSpectrum[binLocal];
              }
            }
          }
          // and keep the mean and rms values as computed
          _lastGoodMean = spectrumSum / nBins;
          _lastGoodRMS = sqrt(spectrumSumSq/nBins - std::pow(_lastGoodMean,2));
        }
      }


      // Now we have a valid spectrum, either the original or
      // replaced; this spectrum is good, so let us do the final
      // bits of post processing reset the spectrumSum, and SumSq,
      // and flatten subtract the bandpass from the data.
      //
      spectrumSum = 0.0;
      spectrumSumSq = 0.0;
      for (unsigned s = 0; s < nSubbands; ++s) 
      {
        long index = t*nChannels*nSubbands+s*nChannels;
        for (unsigned c = 0; c < nChannels; ++c) 
        {
          int binLocal = s*nChannels +c;
          //std::cout << "here1 " << bandPass[binLocal] << " " << dataModel << std::endl;
          // flat bandpass with near zero mean
          I[index+c] -= (bandPass[binLocal] + dataModel);
          //std::cout << "here2" << std::endl;
          spectrumSum += I[index+c];
          //std::cout << "here3" << std::endl;
          spectrumSumSq += I[index+c]*I[index+c];
        }
      }


      // and normalize: bring to zero mean if zerodm is specified or
      // use the running mean if not
      spectrumSum /= nBins; // New meaning of these two variables
      spectrumRMS = sqrt(spectrumSumSq/nBins - std::pow(spectrumSum,2));
      //std::cout << "Values used for normalization" << std::endl;
      //std::cout << spectrumSum << " " << spectrumRMS << std::endl;

      // Avoid nastiness in those first spectra by avoiding divisions
      // by zero
      if (spectrumRMS == 0.0) spectrumRMS = 1.0;
      for (unsigned s = 0; s < nSubbands; ++s) 
        {
          long index = t*nChannels*nSubbands+s*nChannels;
          for (unsigned c = 0; c < nChannels; ++c) 
          {          
            if (_zeroDMing == 1)
            {
              I[index+c] -= _zeroDMing * spectrumSum;
            }
            else
            {
              I[index+c] -= _meanRunAve;
            }
            // it may be better to normalize by the running average RMS,
            // given this is a sensitive operation. For example, an
            // artificially low rms may scale things up
            I[index+c] /= spectrumRMS;
            //std::cout<< "Normalised data: "<< index+c <<" "<<I[index+c] <<std::endl;
            //    I[index+c] /= _rmsRunAve;
            // make sure this division is not introducing signals that
            // you would have clipped
            if (I[index+c] > _crFactor)
            { 
              //std::ofstream f3("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type3.txt", ios::app);
              std::ofstream f3(nameofoutput3, ios::app);
              I[index+c] = 0.0;
              long place = a*averageSamplestoReadFromFile+index + c ;
              f3 << (place) << "\n";
              //SamplesReplacedType3.push_back(place);
              f3.close();
            }
          }
        }
      // The bandpass is flat and the spectrum clean and normalized,
      // so move to the next spectrum and write out some stats:
      unsigned reportStatsEvery = 10 * _maxHistory;
      if (_num == 0) 
      {
        // calculate fractions
        float fractionBadSpectra = 100.0 * (float)_badSpectra / (float)reportStatsEvery;
        float fractionBadChannels = 100.0 * _fractionBadChannels / (float)reportStatsEvery ;

        // if the fraction of bad spectra becomes >99%, then empty the
        // circular buffers and go into learning mode again
        if (fractionBadSpectra > 99.0) 
        {
          _rmsBuffer.resize(0);
          _meanBuffer.resize(0);
          std::cout << "Lost track of the RFI model, retraining.";
        }

        std::cout << std::endl;
        //  std::cout << "# Mean RMS %badSpectra %badChannels : "
        //  << std::endl;
        std::cout <<  "RFIstats: " << _meanRunAve << " " << _rmsRunAve << " "
        << fractionBadSpectra << " " << fractionBadChannels << std::endl;
        std::cout << std::endl;
        std::cout << nBins << std::endl;
        // Reset _bad
        _badSpectra = 0;
        _fractionBadChannels = 0.0;
      }

      // and update the model
      bandPassAve=_meanRunAve;
      bandPassVar=_rmsRunAve;
      ++_num;
      _num = _num % reportStatsEvery;

    } //end of: for (unsigned t = 0; t < nSamples; ++t) loop
    
    outfile.open(outfilename, ios::binary | ios::out | ios::app);
    //outfile.open("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_normalised.fil", ios::binary | ios::out | ios::app);

    outfile.write(reinterpret_cast<char*>(I.data()), I.size()*sizeof(float));
    outfile.close();
  }
 
 tend = time(0); 
 cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;
 std::ifstream f2(nameofoutput2, ios::in);
 fsize = f2.tellg();
 f2.seekg( 0, std::ios::end );
 fsize = f2.tellg() - fsize;
 std::cout << "The size of the file is: " << fsize << std::endl;

 cout << "The number of samples replaces because of type1 RFI is: "<< SamplesReplacedType1.size() <<std::endl;
 // if (SamplesReplacedType1.size()>0)
 // {
 //   std::ofstream f1("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type1.txt", ios::app);
 //   for(vector<float>::const_iterator i = SamplesReplacedType1.begin(); i != SamplesReplacedType1.end(); ++i) 
 //   {
 //      f1 << *i << '\n';
 //   }
 //   f1.close();
 // }
 //cout << "The number of samples replaces because of type2 RFI is: "<< SamplesReplacedType2.size() <<std::endl;
 // if (SamplesReplacedType2.size()>0)
 // {
 //   std::ofstream f2("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type2.txt", ios::app);
 //   for(vector<float>::const_iterator i = SamplesReplacedType2.begin(); i != SamplesReplacedType2.end(); ++i) 
 //   {
 //      f2 << *i << '\n';
 //   }
 //   f2.close();
 // }
 // cout << "The number of samples replaces because of type3 RFI is: "<< SamplesReplacedType3.size() <<std::endl;
 // if (SamplesReplacedType3.size()>0)
 // {
 //   std::ofstream f3("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type3.txt", ios::app);
 //   for(vector<float>::const_iterator i = SamplesReplacedType3.begin(); i != SamplesReplacedType3.end(); ++i) 
 //   {
 //      f3 << *i << '\n';
 //   }
 //   f3.close();
 // }
 return 0;

} //end of: int main(void)



