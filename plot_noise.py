import numpy as np
import matplotlib.pyplot as plt
from rich.pretty import Pretty
from sigpyproc.readers import FilReader
from sys import argv

filname = argv[1]
print(filname)
fil = FilReader(argv[1])
Pretty(fil.header)

end_data = -1
end_data = int(argv[2])

if end_data == -1:
    data = fil.read_block(0,fil.header.nsamples)
else:
    data = fil.read_block(0,end_data)

plt.imshow(data,aspect='auto')
plt.title("Filterbank Data")
plt.show()