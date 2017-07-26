# this is the wrapper function for calling simphy for simulation
# and organize the output gene tree files into one file

# for how to setup the simphy and the control file
# check https://github.com/adamallo/SimPhy/wiki/Manual


simphy.simu<-function(simphy.executable="simphy_lnx64",control.file,remove.individual.gtreefiles=T){
  #run simphy
  simphy.command<-paste('./',simphy.executable, ' -I ',control.file,sep='')
  system(simphy.command)

  #read control file for the output file name

  simphy.contrlfile<-scan(file=control.file,what = 'character',sep = '\n')
  # rl.line<-grep("^-rl ",simphy.contrlfile)
  # rl<-NULL
  #
  # if(length(rl.line)==1){
  #   rl<-as.integer(sub("^-rl \\w+:(\\d+)\\.*$","\\1",simphy.contrlfile[rl.line]))
  # }
  o.line<-grep("^-o ",simphy.contrlfile)
  output.folder=NULL
  if(length(o.line)==1){
    output.folder<-sub("^-o (\\S+) *.*$","\\1",simphy.contrlfile[o.line])
  }

  # if(is.null(rl)){
  #   print("This is the control file:")
  #   print(simphy.contrlfile)
  #   print("Could not find where the number of locus is specified")
  #   print("please enter the number of locus:")
  #   rl=scan(what = "integer",n=1)
  #
  # }
  if(is.null(output.folder)){
    print("This is the control file:")
    print(simphy.contrlfile)
    print("Could not find where the output foldername is specified")
    print("please enter output folder name:")
    output.folder=scan(what = "character",n=1)

  }
  simu.rep.folder<-dir(path = output.folder,pattern = "\\d+")
  if(length(simu.rep.folder)==0 ){
    stop("Could not find the simulated folders")
  }
  for (folder in simu.rep.folder){
    folder.path=paste(output.folder,folder,sep='/')
    clean.simphy.folder(folder.path,remove.individual.gtreefiles)
  }

  return(paste(output.folder,simu.rep.folder,sep='/'))

}
