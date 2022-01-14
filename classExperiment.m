
classdef classExperiment
   properties
      name 
      sample
      rate
      x
      y
      pH
      ramp
      concentration
      loc
      tag
      UFline
      Fline
      FUF
      MF
      MUF
      BF
      BUF
      uy
      ux
      ly
      lx
   end
   methods
       function currExp = getFUF(currExp)
             
           Folded = currExp.MF*currExp.x+currExp.BF;
           Unfolded = currExp.MUF*currExp.x+currExp.BUF;
%             Folded = currExp.Fline;
%             Unfolded =currExp.UFline;
           newFUF = (currExp.y-Folded)./(Unfolded-Folded);

              
           currExp.FUF = newFUF;
           currExp.UFline = Unfolded;
           currExp.Fline = Folded;
       end

       function cellExperiments = autoloadExperiments(temp,filepath) %creates cell array of all 
%             filepath = uigetdir('C:\Users\chris\MATLAB Drive\Summer Project\Data');
            dirs = dir(filepath);
           files = {dirs.name};
           i = 3; %first 2 files are .. and ...
          
           if isequal(files{i},".DS_Store")
               files(i) = [];
           end
            cellExperiments = cell(size(files,2)-2,1);
           while(i <= size(files,2)) %creation of the cell array this loop reads all files in directory
               location = strcat(filepath,'/',files{i}); 
               tempExperiment = load(location);
               tempExperiment = tempExperiment.temp;
               tempExperiment.loc = strcat(filepath,'/');
               
               
               
               cellExperiments{i-2} = tempExperiment;
               
               i = i+1;
           end
       end
       



   end
end