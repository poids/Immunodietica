#Script checks for required packages and installs them if they are not found
#By Vasco Morais

packages=c('dplyr', 'ggplot2', 'stringr') #Change which packages are required here

install_required_packages=function(packages) {
  #Checks which packages are already installed
  installed_packages=(packages %in% row.names(installed.packages()))
  
  if (FALSE %in% installed_packages) {
    install_packages=packages[!installed_packages]
    print(paste('Installing', install_packages))
    install.packages(install_packages, dependencies = TRUE, repos = 'http://cran.us.r-project.org')
  } else {
    print('All required R libraries/packages found!')
  }
}

install_required_packages(packages)