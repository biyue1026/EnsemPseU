# EnsemPseU
  
EnsemPseU is an ensemble approach to identify pseudouridine sites. It can be implemented by running the ‘EnsemPseU_pred.py’. 
The parameters of EnsemPseU_pred.py are:
* **--species**: the species that you want to predict.You can choose 'H','M' and 'S', which represent *H. sapiens*, *M. musculus* and                          *S. cerevisiae*, respectively.
  
* **--input**: the input RNA sequence file with FASTA format. Please note that the length of input sequence must be 21-bp when the 
         species is *H. sapiens* or *M. musculus*, and 31-bp when the species is *S. cerevisiae*.
           
* **--output**: the output file name. The prediction results would be saved in csv format.
