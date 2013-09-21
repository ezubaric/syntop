=============================================
to generate toy data or parse wackypedia data
=============================================

1. set the flag GENERATE in parse_data.sh to true
2. if parsing wackpedia data, set proper VOCAB_SIZE, DOC_LIMIT
3. if generating toy data, set the parameter for python script create_synthetic.py
   for example, python scripts/create_synthetic.py toy/toy 20 10
   toy/toy is the desired prefix of the data
   20 is the document limit
   10 is the average number of words per document
4. do make sure there are no empty docuemnt/tokens in the .doc files, this would cause a parsing error.

Attention: if you are running the code in local mode, you probably need to set
model_name to proper value in syntop_mapper.cpp and syntop_reducer.cpp, i.e.,
params->set_model_name(DIRECTORY)

==============================
to run syntop in a local mode.
==============================

1. set the flag LOCAL in command.sh to true
2. select the dataset you want to use, and then select the parameters.
   DATASET="wackypedia"
   VOCAB_SIZE=5050
   DOC_LIMIT=983

   DATASET="toy"
   VOCAB_SIZE=45
   DOC_LIMIT=20

   DATASET="test_data" (this is probably the best synthetic data we could work on)
   VOCAB_SIZE=45
   DOC_LIMIT=500
3. open cli and type in command "bash command.sh"
4. syntop output will be saved to output*.txt
5. model will be saved to DATASET_param folder
