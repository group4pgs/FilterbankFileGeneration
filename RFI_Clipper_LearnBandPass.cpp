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
#include <sstream>
#include <fftw3.h>
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

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::string getexepath()
{
    char *path = NULL;
    path = getcwd(NULL, 0); // or _getcwd
    if ( path != NULL)
    return path;
}


/**
*********************************************************************
2. Beginning of the main function.
*********************************************************************
**/
int main(int argc, char * argv [])
{
   time_t tstart, tend, tstartDedisperse, tendDedisperse , tstartAutocorr, tendAutocorr; 
   tstart = time(0);

  /**
  *********************************************************************
  3. Set the default values before processing
  *********************************************************************
  **/
  int nSubbands   = 1;
  unsigned goodSamples = 0;
  float dataChunksInSec;
  float spectrumRMS = 0.0;
  float spectrumSum = 0.0;
  char  fileNumber;
  std::string outfilename, directory, infilename;
  float _crFactor, _srFactor, averageFilterLength;
  float _rmsRunAve  = 0.0;
  float _meanRunAve = 0.0;
  float _useMeanOverRMS = 0;
  float _zeroDMing = 1;
  float k = 4; // degrees of freedom
  int nbands;
  float meanMinusMinimum = k / sqrt(2.0*k);
  std::vector<float> _lastGoodSpectrum, _dataBandPass;
  int _badSpectra(0), _badChannels(0);
  int _current(0); // history pointer
  int NumberOfFlaggedSamples = 0;
  //std::pair<int, int> apair;
  //vector<std::pair<int, int>> SamplesReplacedType1, SamplesReplacedType2, SamplesReplacedType3;
  boost::circular_buffer<float> _meanBuffer, _rmsBuffer;
  std::vector<unsigned char> n(1000);
  std::vector<unsigned char> tsampRaw(8);
  std::vector<unsigned char> nbitsRaw(4);
  std::vector<unsigned char> nchansRaw(4);
  double tsamp,_meanOverRMS;
  float _fractionBadChannels(0.0);
  int nChannels, nbits, headerSize, plek, _remainingZeros;
  bool _num = 0;

  /**
  *********************************************************************
  3. Parse the input arguments.
  *********************************************************************
  **/

  if (cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "-help")) 
  {
      std::cout<< "Options for the RFI clipper are:"<< std::endl;
      std::cout<< "\n";
      std::cout<< "-inputFile [-]\t\t" << "The name of the filterbank file on which the RFI Clipper will run." <<std::endl;
      std::cout<< "-outputFile [-]\t\t"<< "The name of the output file. (default: 'output.fil')" <<std::endl;
      std::cout<< "-outputDirectory [-]\t"<< "Specify the output directory if different to the input directory." <<std::endl;
      std::cout<< "-fileNumber [-]\t\t"<< "Specify the file number in case they have the same name but different numbers." <<std::endl;
      std::cout<< "-filterLength [-]\t"<< "Specify the moving average length in seconds (def 1 sec)." <<std::endl;
      std::cout<< "-dataChunksInSec [-]\t"<< "Process the whole data file in smaller chunks. The size of the smaller chunks should be specified in seconds (def 1 sec)." <<std::endl;
      std::cout<< "-crFactor [-]\t\t"<< "Margin of tolerance for bright channels (def 10)." <<std::endl;
      std::cout<< "-srFactor [-]\t\t"<< "Margin of tolerance for bright spectrums (all the channels for a specific time instance) (def 4)." <<std::endl;
      std::cout<< "-subDivisions [-]\t"<< "Specify the number of subdivisions of the bandpass for building up statistics necessary to determine clipping." <<std::endl;
      std::cout<< "\t\t\tThis is to minimize the computational overhead of the overall function (default number of subdivisions = 8)." <<std::endl;
      std::cout<< "\n";
      std::exit(0);

  }

  if (cmdOptionExists(argv, argv+argc, "-crFactor"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-crFactor");
    _crFactor = atof(temp);
  }

  else if (cmdOptionExists(argv, argv+argc, "-crfactor"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-crfactor");
    _crFactor = atof(temp);
  }
  else
  {
    _crFactor = 10; // channel rejection RMS 
  }

  if (cmdOptionExists(argv, argv+argc, "-srFactor"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-srFactor");
    _srFactor = atof(temp);
  }

  else if (cmdOptionExists(argv, argv+argc, "-srfactor"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-srfactor");
    _srFactor = atof(temp);
  }
  else
  {
    _srFactor = 4; // spectrum rejection RMS
  }


  if (cmdOptionExists(argv, argv+argc, "-inputFile"))
  {
      infilename = getCmdOption(argv, argv + argc, "-inputFile");
  }
  else if (cmdOptionExists(argv, argv+argc, "-inputfile"))
  {
      infilename = getCmdOption(argv, argv + argc, "-inputfile");
  }
  else
  {
       cout << "No input file was specified." << std::endl;
       exit(0);
  }

  if (cmdOptionExists(argv, argv+argc, "-outputFile"))
  {
    outfilename = getCmdOption(argv, argv + argc, "-outputFile");
  }

  else if (cmdOptionExists(argv, argv+argc, "-outputfile"))
  {
    outfilename = getCmdOption(argv, argv + argc, "-outputfile");
  }
  else
  {
    cout << "No output file was specified. The default name will be used: 'output.fil'" << std::endl;
    outfilename = "output.fil";
  }

  if (cmdOptionExists(argv, argv+argc, "-outputDirectory"))
  {
    directory = getCmdOption(argv, argv + argc, "-outputDirectory");
  }

  else if (cmdOptionExists(argv, argv+argc, "-outputdirectory"))
  {
    directory = getCmdOption(argv, argv + argc, "-outputdirectory");
  }
  else
  {
      cout << "No output directory was specified." << std::endl;
      directory = getexepath();
  }
  

  if (cmdOptionExists(argv, argv+argc, "-fileNumber"))
  {
    fileNumber = *getCmdOption(argv, argv + argc, "-fileNumber");
  }

  else if (cmdOptionExists(argv, argv+argc, "-filenumber"))
  {
    fileNumber = *getCmdOption(argv, argv + argc, "-filenumber");
  }
  else
  {
    cout << "No file number was specified. The default value of '0' will be used." << std::endl;
    fileNumber = '0';
  }

  if (cmdOptionExists(argv, argv+argc, "-filterLength"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-filterLength");
    averageFilterLength = atof(temp);
  }

  else if (cmdOptionExists(argv, argv+argc, "-filterlength"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-filterlength");
    averageFilterLength = atof(temp);
  }
  else
  {
    cout << "No file number was specified. The default value of '1' will be used." << std::endl;
    averageFilterLength = 1.0;
  }

  if (cmdOptionExists(argv, argv+argc, "-dataChunksInSec"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-dataChunksInSec");
    dataChunksInSec = atof(temp);
  }

  else if (cmdOptionExists(argv, argv+argc, "-datachunksinsec"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-datachunksinsec");
    dataChunksInSec = atof(temp);
  }
  else
  {
    dataChunksInSec = 5; // channel rejection RMS 
  }

  if (cmdOptionExists(argv, argv+argc, "-subDivisions"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-subDivisions");
    nbands = atoi(temp);
  }

  else if (cmdOptionExists(argv, argv+argc, "-subdivisions"))
  {
    char * temp = getCmdOption(argv, argv + argc, "-subdivisions");
    nbands = atoi(temp);
  }
  else
  {
    nbands = 8; // channel rejection RMS 
  }

  /**
  ****************************************************************************
  4. Read the header of the SIGPROC file and write the header to the new file
  ****************************************************************************
  **/

  // Determine the size of the whole file
  std::streampos fsize = 0;
  ifstream infile;
  cout << infilename << endl;
  //infile.open("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0.fil", ios::binary | ios::in);
  infile.open(infilename, ios::binary | ios::in);
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
  std::cout<< "Number of channels are: "<< nBins <<"\n";

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

  float _maxHistory  = averageFilterLength/tsamp;  // Filter length 
  _meanBuffer.set_capacity(_maxHistory);
  _rmsBuffer.set_capacity(_maxHistory);
  int tobs = (float(fsize)-headerSize)/(nChannels*nSubbands)*tsamp;
  int nSamples = (float(fsize)-headerSize)/(nbits/8.0)/(nChannels*nSubbands)-1;
  int NumberOfChunks  = floor(tobs/dataChunksInSec);
  int NumberOfSamplesPerChunk = floor(dataChunksInSec/tsamp);
  int averageSamplestoReadFromFile = (nChannels*nSubbands*NumberOfSamplesPerChunk);


  // std::string nameofoutput1 =  directory + std::string(fileNumber) + std::string("_RFI_Type1.txt"); 
  // std::string nameofoutput2 =  directory + std::string(fileNumber) + std::string("_RFI_Type2.txt");
  // std::string nameofoutput3 =  directory + std::string(fileNumber) + std::string("_RFI_Type3.txt"); 

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
      _dataBandPass.resize(nBins,0.0);
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
      float spectrumSum = 0.0;
      bool goodSpectrum = true;
      float spectrumSumSq = 0.0;
      // Split the spectrum into n bands for the purpose of matching
      // the spectrum to the model
      std::vector<float> goodChannels(nbands,0.0);
      unsigned channelsPerBand = nBins / nbands;
      // find the minima of I in each band, and the minima of the
      // model and compare
      std::vector<float> miniData(nbands,1e6);
      std::vector<float> miniModel(nbands,1e6);
      std::vector<float> dataMinusModel(nbands);
      std::vector<float> bandSigma(nbands);
      std::vector<float> bandMean(nbands,0.0);
      std::vector<float> bandMeanSquare(nbands,0.0);


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

            //if (I[index+c] < miniData[band]) miniData[band] = I[index+c];
            //if (bandPass[binLocal] < miniModel[band]) miniModel[band] = bandPass[binLocal];
            bandMean[band] += I[index+c];
            bandMeanSquare[band] += (I[index+c]*I[index+c]);
          }
        }


        // Now find the distances between data and model and the RMSs in
        // each band

        for (unsigned b = 0; b < nbands; ++b)
        {
          //dataMinusModel[b] = miniData[b] - miniModel[b];
          bandMean[b] /= channelsPerBand; 
          bandSigma[b] = sqrt(bandMeanSquare[b]/channelsPerBand - std::pow(bandMean[b],2));
          //std::cout << bandSigma[b] << " " << bandMean[b] << " " <<
          //bandMeanSquare[b] << std::endl;
        }


		// Recompute _dataBandPass as a sequence of flat segments,
		// based on the average total power per band, also averaged
		// over the training samples
        for (unsigned s = 0; s < nSubbands; ++s) 
        {
          long index = t*nChannels*nSubbands+s*nChannels;
          for (unsigned c = 0; c < nChannels ; ++c) 
          {
            int binLocal = s*nChannels +c;
            unsigned band = (int)binLocal / channelsPerBand;
	    	_dataBandPass[binLocal] = (_dataBandPass[binLocal] * _rmsBuffer.size() + bandMean[band])/(_rmsBuffer.size()+1);
          }
        }

        // Assume the minimum bandSigma to be the best estimate of this
        // spectrum RMS
        spectrumRMS = *std::min_element(bandSigma.begin(), bandSigma.end());
        
        // Take the median of dataMinusModel to determine the distance from the model
        //std::nth_element(dataMinusModel.begin(), dataMinusModel.begin()+dataMinusModel.size()/2, dataMinusModel.end());
        //dataModel = (float)*(dataMinusModel.begin()+dataMinusModel.size()/2);
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
          //dataModel = 0.0;
        }

        // and update the model
        //dataModel = dataModel + meanMinusMinimum * _rmsRunAve;
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
        // dataModel = 0.0;
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
          if (I[index+c] - _dataBandPass[binLocal] > margin) 
          {
            // clipping this channel to values from the last good
            // spectrum
            /*std::cout << "clipping channel " << I[index+c]
              << " " << dataModel << " " <<  bandPass[binLocal]
              << " " << margin << std::endl;*/
            //The following is for polarization
            I[index + c] = _lastGoodSpectrum[binLocal];
            // apair.first = a*averageSamplestoReadFromFile+index;
            // apair.second = c;
            // SamplesReplacedType1.push_back(apair);
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
      float tempMean = std::accumulate(_dataBandPass.begin(), _dataBandPass.end(), 0.0)/nBins;

      // Check if more than 20% of the channels in each band were
      // bad. If so in more than half of the bands, 4 in this case,
      // keep record. Also, if one band is completely gone, or less
      // than 80% of the total survive, keep record.
      unsigned badBands = 0;
      for (unsigned b = 0; b < nbands; ++b)
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
      float spectrumRMStolerance = _srFactor * _rmsRunAve/sqrt(nBins); 

      //Now check, if spectrumSum - model > tolerance, declare this
      //time sample useless, replace its data and take care of the
      //running averages, also cut spectra where badBands >= 4, see above

      if (spectrumSum - _meanRunAve > spectrumRMStolerance || badBands >= 4) 
      {
        // we need to remove this entire spectrum
        for (unsigned s = 0; s < nSubbands; ++s) 
        {
            long index = t*nChannels*nSubbands+s*nChannels;
            for (unsigned c = 0; c < nChannels; ++c) 
            {
              int binLocal = s*nChannels +c;
              I[index + c] = _lastGoodSpectrum[binLocal];
              // apair.first = a*averageSamplestoReadFromFile+index;
              // apair.second = c;
              // SamplesReplacedType2.push_back(apair);
            }
        }
        // goodSpectrum = false;  
        // keep a record of bad spectra
        ++_badSpectra;

        // now remove the last samples from the running average
        // buffers and replace them with the last good values.
        //
        _meanBuffer.pop_back();
        _rmsBuffer.pop_back();

        _meanBuffer.push_back(_meanBuffer.back());
        _rmsBuffer.push_back(_rmsBuffer.back());
      }
      // else keep a copy of the original spectrum, as it is good
      else
      {
    		if (_current < 10000)
        // if (_current < _maxHistory)
    		{
    		    float dataBandPassMean;
    		    float offset;
    		    ++_current;

    		    for (unsigned s = 0; s < nSubbands; ++s) 
            	{
                	long index = t*nChannels*nSubbands+s*nChannels;
                	for (unsigned c = 0; c < nChannels; ++c) 
                	{
                  		int binLocal = s*nChannels +c;
    					_dataBandPass[binLocal] = (_dataBandPass[binLocal] * (_current - 1) + I[index+c]) / _current;
                	}
            	}
    		    // Now place the band pass at the right level,
    		    // i.e. subtract its mean and add the running average of
    		    // the data
    		    dataBandPassMean = std::accumulate(_dataBandPass.begin(), _dataBandPass.end(), 0.0)/nBins;
    		    offset = _meanRunAve - dataBandPassMean;
    		    transform(_dataBandPass.begin(), _dataBandPass.end(), _dataBandPass.begin(), bind2nd(std::plus<float>(), offset));
    		    
    		    // on the last iteration, write out the bandpass and keep the mean
    		    if (_current == _maxHistory - 1)
    		    {
    				// Copy the level of the data bandpass into a global
    				// the following two lines are just a test
    				/*
    				std::cout << "BandPass formed "<< meanBand << 
    				  " " << _meanRunAve << std::endl;
    				*/
    				//_dataBandPassLevel = _meanRunAve;
    				// Write out the bandpass to a file
    				std::stringstream BandPassOutputName;
    				BandPassOutputName << "bandpass_" << nSubbands << ".dat";
    				std::ofstream output_file(BandPassOutputName.str().c_str());
    				std::ostream_iterator<float> output_iterator(output_file, "\n");
    				std::copy(_dataBandPass.begin(), _dataBandPass.end(), output_iterator);
    		    }
    		}
    		else
    		{
    		    float offset, dataBandPassMean;
    		    // Now place the band pass at the right level,
    		    // i.e. subtract its mean and add the running average of
    		    // the data
    		    dataBandPassMean = std::accumulate(_dataBandPass.begin(), _dataBandPass.end(), 0.0)/nBins;
    		    //	    offset = _meanRunAve - _dataBandPassLevel;
    		    offset = _meanRunAve - dataBandPassMean;
    		    transform(_dataBandPass.begin(), _dataBandPass.end(), _dataBandPass.begin(), bind2nd(std::plus<float>(), offset));
    		}
		// At this stage I have a bandpass from the data: TEST

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
          //I[index+c] -= (bandPass[binLocal] + dataModel);
          I[index+c] -= _dataBandPass[binLocal]; 
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
              I[index+c] = 0.0;
              // apair.first = a*averageSamplestoReadFromFile+index;
              // apair.second = c;
              // SamplesReplacedType3.push_back(apair);
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
          _current = 0;
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
      ++_num;
      _num = _num % reportStatsEvery;

    } //end of: for (unsigned t = 0; t < nSamples; ++t) loop

    outfile.open(outfilename, ios::binary | ios::out | ios::app);
    outfile.write(reinterpret_cast<char*>(I.data()), I.size()*sizeof(float));
    outfile.close();
  }
 tend = time(0); 
 cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

 
 // if (SamplesReplacedType1.size()>0)
 // {
 //   cout << "The number of samples replaces because of type1 RFI is: "<< SamplesReplacedType1.size() <<std::endl;
 //   std::ofstream f1("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type1.txt", ios::app);
 //   for(vector<pair<int, int>>::iterator it = SamplesReplacedType1.begin(); it != SamplesReplacedType1.end(); ++it) 
 //   {
 //   	  apair = *it;
 //      f1 << apair.first << "," << apair.second << '\n';
 //   }
 //   f1.close();
 // }
 // if (SamplesReplacedType2.size()>0)
 // {
 //   cout << "The number of samples replaces because of type2 RFI is: "<< SamplesReplacedType2.size() <<std::endl;
 //   std::ofstream f2("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type2.txt", ios::app);
 //   for(vector<pair<int,int>>::iterator it = SamplesReplacedType2.begin(); it != SamplesReplacedType2.end(); ++it) 
 //   {
 //   	  apair = *it;
 //      f2 << apair.first << "," << apair.second << '\n';
 //   }
 //   f2.close();
 // }
 // if (SamplesReplacedType3.size()>0)
 // {
 //    cout << "The number of samples replaces because of type3 RFI is: "<< SamplesReplacedType3.size() <<std::endl;
 //   std::ofstream f3("/home/elmarie/Journal Publications/RFI Documentation/Data/NonStat/SNInject0_RFI_Type3.txt", ios::app);
 //   for(vector<pair<int,int>>::iterator it = SamplesReplacedType3.begin(); it != SamplesReplacedType3.end(); ++it) 
 //   {
 //   	  apair = *it;
 //      f3 << apair.first << "," << apair.second << '\n';
 //   }
 //   f3.close();
 // }

 return 0;
} //end of: int main(void)



