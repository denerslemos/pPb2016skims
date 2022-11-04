universe   = vanilla
getenv     = True
executable = sub_skim.sh
+JobFlavour           = "tomorrow"
requirements = (OpSysAndVer =?= "CentOS7")
RequestCpus = 2
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"

log        = cond/test.pPbpthat15.log
output     = cond/test.pPbpthat15.out
error      = cond/test.pPbpthat15.err
arguments = listoffiles_pPb_MC_pthat15_pgoing.txt HiForestAODskim_pthat15_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat30.log
output     = cond/test.pPbpthat30.out
error      = cond/test.pPbpthat30.err
arguments = listoffiles_pPb_MC_pthat30_pgoing.txt HiForestAODskim_pthat30_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat50.log
output     = cond/test.pPbpthat50.out
error      = cond/test.pPbpthat50.err
arguments = listoffiles_pPb_MC_pthat50_pgoing.txt HiForestAODskim_pthat50_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat80.log
output     = cond/test.pPbpthat80.out
error      = cond/test.pPbpthat80.err
arguments = listoffiles_pPb_MC_pthat80_pgoing.txt HiForestAODskim_pthat80_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat120.log
output     = cond/test.pPbpthat120.out
error      = cond/test.pPbpthat120.err
arguments = listoffiles_pPb_MC_pthat120_pgoing.txt HiForestAODskim_pthat120_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat170.log
output     = cond/test.pPbpthat170.out
error      = cond/test.pPbpthat170.err
arguments = listoffiles_pPb_MC_pthat170_pgoing.txt HiForestAODskim_pthat170_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat220.log
output     = cond/test.pPbpthat220.out
error      = cond/test.pPbpthat220.err
arguments = listoffiles_pPb_MC_pthat220_pgoing.txt HiForestAODskim_pthat220_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat280.log
output     = cond/test.pPbpthat280.out
error      = cond/test.pPbpthat280.err
arguments = listoffiles_pPb_MC_pthat280_pgoing.txt HiForestAODskim_pthat280_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat370.log
output     = cond/test.pPbpthat370.out
error      = cond/test.pPbpthat370.err
arguments = listoffiles_pPb_MC_pthat370_pgoing.txt HiForestAODskim_pthat370_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat460.log
output     = cond/test.pPbpthat460.out
error      = cond/test.pPbpthat460.err
arguments = listoffiles_pPb_MC_pthat460_pgoing.txt HiForestAODskim_pthat460_pgoing_MC.root true 
queue

log        = cond/test.pPbpthat540.log
output     = cond/test.pPbpthat540.out
error      = cond/test.pPbpthat540.err
arguments = listoffiles_pPb_MC_pthat540_pgoing.txt HiForestAODskim_pthat540_pgoing_MC.root true 
queue


log        = cond/test.Pbppthat15.log
output     = cond/test.Pbppthat15.out
error      = cond/test.Pbppthat15.err
arguments = listoffiles_pPb_MC_pthat15_Pbgoing.txt HiForestAODskim_pthat15_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat30.log
output     = cond/test.Pbppthat30.out
error      = cond/test.Pbppthat30.err
arguments = listoffiles_pPb_MC_pthat30_Pbgoing.txt HiForestAODskim_pthat30_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat50.log
output     = cond/test.Pbppthat50.out
error      = cond/test.Pbppthat50.err
arguments = listoffiles_pPb_MC_pthat50_Pbgoing.txt HiForestAODskim_pthat50_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat80.log
output     = cond/test.Pbppthat80.out
error      = cond/test.Pbppthat80.err
arguments = listoffiles_pPb_MC_pthat80_Pbgoing.txt HiForestAODskim_pthat80_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat120.log
output     = cond/test.Pbppthat120.out
error      = cond/test.Pbppthat120.err
arguments = listoffiles_pPb_MC_pthat120_Pbgoing.txt HiForestAODskim_pthat120_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat170.log
output     = cond/test.Pbppthat170.out
error      = cond/test.Pbppthat170.err
arguments = listoffiles_pPb_MC_pthat170_Pbgoing.txt HiForestAODskim_pthat170_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat220.log
output     = cond/test.Pbppthat220.out
error      = cond/test.Pbppthat220.err
arguments = listoffiles_pPb_MC_pthat220_Pbgoing.txt HiForestAODskim_pthat220_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat280.log
output     = cond/test.Pbppthat280.out
error      = cond/test.Pbppthat280.err
arguments = listoffiles_pPb_MC_pthat280_Pbgoing.txt HiForestAODskim_pthat280_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat370.log
output     = cond/test.Pbppthat370.out
error      = cond/test.Pbppthat370.err
arguments = listoffiles_pPb_MC_pthat370_Pbgoing.txt HiForestAODskim_pthat370_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat460.log
output     = cond/test.Pbppthat460.out
error      = cond/test.Pbppthat460.err
arguments = listoffiles_pPb_MC_pthat460_Pbgoing.txt HiForestAODskim_pthat460_Pbgoing_MC.root true 
queue

log        = cond/test.Pbppthat540.log
output     = cond/test.Pbppthat540.out
error      = cond/test.Pbppthat540.err
arguments = listoffiles_pPb_MC_pthat540_Pbgoing.txt HiForestAODskim_pthat540_Pbgoing_MC.root true 
queue