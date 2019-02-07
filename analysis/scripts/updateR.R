#update R

# installing/loading the package:
if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr

#Using the package:
#This will start the updating process of your R installation.  
#It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.
updateR() 

#https://stackoverflow.com/questions/13656699/update-r-using-rstudio