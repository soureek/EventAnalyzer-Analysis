EventAnalyzer for Analysis with CMS Open Data

Follow [CMS Open Data Workshop](https://cms-opendata-workshop.github.io/workshopwhepp-lesson-docker/03-docker-for-cms-opendata/index.html) to set up the 
docker container.

**Setup instrauctions:**
```
cd CMSSW_7_6_7/src/
cmsenv
git clone git@github.com:soureek/EventAnalyzer-Analysis.git EventAnalyzer/Analysis
scram b -j 9
```

**Running instructions:**
  
  - `test/runAnalysis_cfg.py` is the configuration file to do `cmsRun` and it creates an output file `OutFile_MC.root` or `OutFile_data.root`. The Analysis code is defined in `plugins/Analysis.cc`
  
  - Set recommended GlobalTag for MC and Data separately using `options.globalTag` in the python configuration or with command-line arguments with `cmsRun`
    as shown below.
  
  - For MC:   
    ```
    cd $CMSSW_BASE/src/EventAnalyzer/Analysis/test
    cmsRun runAnalysis_cfg.py isMC=True maxEvts=1000
    ```
  
  - For Data:
    ```
    cd $CMSSW_BASE/src/EventAnalyzer/Analysis/test
    cmsRun runAnalysis_cfg.py isMC=False maxEvts=1000
    ```
