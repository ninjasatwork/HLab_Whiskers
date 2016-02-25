function loadSUData(cellNum, SU);
 dirName = uigetdir('', 'Select base directory')
 SU_ConDir = ([dirName filesep 'ConTA' filesep]);
 SU_TDir =   ([dirName filesep 'TrialArrays' filesep]);

 display(['Loading '  SU.trialArrayName{cellNum}])

load([SU_TDir SU.trialArrayName{cellNum}]);
load([SU_ConDir SU.contactsArrayName{cellNum}]);
assignin('base','T',T)
assignin('base','contacts',contacts)
assignin('base','params',params)




