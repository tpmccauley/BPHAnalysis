# BPHAnalysis

## Running in docker

If you don't have docker installed already, instructions are [here](https://docs.docker.com/install/).

For a more recent CMSSW release like CMSSW_9_4_4, run a docker image:

`
docker run --name bph -it clelange/cmssw:9_4_4 /bin/bash
`

where we use the already-made CMSSW_9_4_4 release from [dockerhub](https://hub.docker.com/r/clelange/cmssw/tags)

Note: as described [here](https://github.com/clelange/cmssw-docker/) this will install a stand-alone CMSSW image so it will be
at least a few GBs and may take a while (but only has to be done once).

Once done, you should see the commmand prompt (`docker run -i` means run interactively):

`
cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $
`

Then clone this repository :

```
cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $ git clone https://github.com/tpmccauley/BPHAnalysis.git
Cloning into 'BPHAnalysis'...
remote: Enumerating objects: 45, done.
remote: Counting objects: 100% (45/45), done.
remote: Compressing objects: 100% (34/34), done.
remote: Total 1327 (delta 15), reused 32 (delta 11), pack-reused 1282
Receiving objects: 100% (1327/1327), 324.45 KiB | 0 bytes/s, done.
Resolving deltas: 100% (791/791), done.
[16:40:29] cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $ ls
BPHAnalysis
```

And build:

```
cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $ scram b
Reading cached build data
>> Local Products Rules ..... started
>> Local Products Rules ..... done
>> Building CMSSW version CMSSW_9_4_4 ----
>> Entering Package BPHAnalysis/uty
>> Leaving Package BPHAnalysis/uty
>> Package BPHAnalysis/uty built
>> Entering Package BPHAnalysis/SpecificDecay
>> Creating project symlinks
Entering library rule at src/BPHAnalysis/SpecificDecay/plugins
>> Compiling edm plugin /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/SpecificDecay/plugins/BPHWriteSpecificDecay.cc 
>> Compiling edm plugin /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/SpecificDecay/plugins/BPHHistoSpecificDecay.cc 
>> Entering Package BPHAnalysis/RecoDecay
Entering library rule at BPHAnalysis/RecoDecay
>> Compiling  /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/RecoDecay/src/BPHRecoCandidate.cc 
>> Compiling  /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/RecoDecay/src/BPHKinematicFit.cc 
>> Compiling  /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/RecoDecay/src/BPHPlusMinusCandidate.cc 
>> Compiling  /home/cmsbld/CMSSW_9_4_4/src/BPHAnalysis/RecoDecay/src/BPHRecoSelect.cc 

...
```

You can exit the container and get back to your original command prompt:

`
cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $ exit
`

To check the status of the container:

`
docker ps -a
`

```
CONTAINER ID        IMAGE                  COMMAND                  CREATED             STATUS                          PORTS               NAMES
3eebec38509d        clelange/cmssw:9_4_4   "/bin/sh -c /bin/zsh…"   35 minutes ago      Exited (0) About a minute ago                       bph
```

To start the container again and continue your work:

`
docker start -i bph
`

which will restart your container

```
cmsbld@3eebec38509d ~/CMSSW_9_4_4/src $ ls
BPHAnalysis
```


### Grid proxy

If you want to get a grid proxy you can mount your `.globus` directory by add `-v ~/.globus:/home/cmsusr/.globus` to your `docker run` command:

`
docker run -v ~/.globus:/home/cmsusr/.globus --name bph -it clelange/cmssw:9_4_4 /bin/bash
`

Then in the running container you can run 

```
cmsbld@a5f12b5b3b6a ~/CMSSW_9_4_4/src $ voms-proxy-init -key /home/cmsusr/.globus/userkey.pem -cert /home/cmsusr/.globus/usercert.pem --voms cms 
```

### Running examples

There are example configuration files in 

```
BPHAnalysis/SpecificDecay/test
BPHAnalysis/RecoDecay/test
```

`BPHAnalysis/SpecficDecay/test/cfg_recoWahis.py` takes as its input AOD files and outputs a skim of selected events in `reco.root`, histograms in `his.root`, and summary text files. 

### Copy files from docker container to host

`
docker cp bph:/home/cmsbld/CMSSW_9_4_4/src/his.root
`

### Running GUIs 

#### Mac OSX

Follow the instructions [here](https://cntnr.io/running-guis-with-docker-on-mac-os-x-a14df6a76efc) which can be boiled down to installing `socat` and `xquartz`:

```
brew install socat
brew install xquartz
```

Run `socat` and open `Xquartz`:

```
socat TCP-LISTEN:6000,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\"
open -a Xquartz
```

Fetch the ip of the host:

```
ifconfig en0
``` 

and note the ip address of `inet`.

Then specifiy the ip on the command line as `DISPLAY` when running:

```
docker run -e DISPLAY=194.12.179.202:0 --rm -it clelange/cmssw:9_4_4 /bin/bash
```

