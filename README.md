# pPb skim code using HTCondor

Code to produce jets, tracks and event plane skims from the CMS HiForest

## Intructions

Setup CMSSW (just for root versioning)
```
export SCRAM_ARCH=slc7_amd64_gcc10
cmsrel CMSSW_12_5_0
cd CMSSW_12_5_0/src
cmsenv
```
Inside of the src folder, download the code using
```
git clone git@github.com:denerslemos/pPbskims.git
cd pPbskims
mkdir cond
```
Before compile the code you must check the [sub_skim.sh](https://github.com/denerslemos/pPbskims/blob/main/sub_skim.sh) lines 4 (CMSSW/src) and 6 (.../pPbskims) and replace by your own folders. You also must replace [line 30 of pPbSkim.C](https://github.com/denerslemos/pPbskims/blob/main/pPbSkim.C#L30) to your own EOS path.

Once this steps are done you can compile the code with
```
g++ -O2 pPbSkim.C `root-config --libs` `root-config --cflags` -o pPbSkim
```
This will create the executable: ```pPbSkim``` 

After that you will need your VOMS certificate, do it using
```
voms-proxy-init -rfc -voms cms --out voms_proxy.txt --hours 200
```
that creates a certificate file valid for 200 hours: ```voms_proxy.txt```

Now you can submit the condor jobs using the python script, [```HTCondor_submit.py```](https://github.com/denerslemos/pPbskims/blob/main/HTCondor_submit.py):

```
python HTCondor_submit.py -i input_text_file -o output_name_file -m X -n Y -s Z
```

- input_text_file: is the text file (use it without the .txt extension) with inputs and can be found in the folders [MC_PYTHIA_SAMPLES](https://github.com/denerslemos/pPbskims/tree/main/MC_PYTHIA_SAMPLES) or [DATA_SAMPLES](https://github.com/denerslemos/pPbskims/tree/main/DATA_SAMPLES) each .root input will be a job

- output_name_file: output file name (use it without the .root extension), it will automatically include a counter for each input

- X: 0 for data and 1 for MC

- Y: 0 for no multiplicity cut (mostly MC or jet samples), 1 for MB [0,185], 2 for HM185 [185,250] and 3 for HM250 [250,inf]

- Z: name for the submission files, I have used HTcondor_sub_ + some information from the sample, pthat, MB, ... + pgoing or Pbgoing.

It will automatically include a counter for each input
