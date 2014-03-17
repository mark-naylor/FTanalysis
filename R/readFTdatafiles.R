
readVermishTestFile <- function(dir="/Users/mnaylor/", filename="FTdata_Ver.txt", ending="", seperator=","){
  
  
  tmp <- file(paste(dir,filename, ending ,sep=""), "rt")   ;  details = readLines(tmp, 1)  ;  close(tmp)
  header <- read.table( paste(dir,filename, ending, sep=""), skip=1, sep= seperator, nrows=2)
  data <- read.table( paste(dir,filename, ending, sep=""), skip=3, sep= seperator)
  
  print("Using Zeta method")
  # Use Zeta method
  Zeta = header[1,1]				# [yr cm^2]
  StErrZeta = header[1,2]			# [yr cm^2]
  REZeta = StErrZeta / Zeta		# [%]
  SqSize = NULL 			# [cm^2]
  
  RhoD = abs (header[2,1])		# [tracks/cm^2]
  StErrRhoD = header[2,2]
  RERhoD = StErrRhoD/RhoD		# [%]
  
  # EffectiveUContent = header[1,3]	# [ppm]
  
  Zeta0 = Zeta
  REZeta0 = REZeta
  RhoD0 = RhoD
  RERhoD0 = RERhoD
  
  K = (Zeta0 * RhoD0) / (Zeta * RhoD)
  # SqNum is known as Dpar in Vers radialplotter
  inputFileData = list(file=filename, details=details, nS = data[,1] , nI = data[,2], SqNum = data[,3], RhoD = RhoD, Zeta  = Zeta, REZeta = REZeta , RERhoD = RERhoD , SqSize = SqSize , StErrZeta = StErrZeta, K=K, StErrRhoD= StErrRhoD)
  
  return(inputFileData)
}



readTarimW2File <- function(dir="/Users/mnaylor/", filename="TarimFTW2.csv", ending="", seperator=","){
  
  tmp <- file(paste(dir,filename, ending ,sep=""), "rt")   ;  details = readLines(tmp, 1)  ;  close(tmp)
  header <- read.table( paste(dir,filename, ending, sep=""), skip=1, sep= seperator, nrows=2)
  data <- read.table( paste(dir,filename, ending, sep=""), skip=3, sep= seperator)
  
  print("Using Zeta method")
  # Use Zeta method
  Zeta = header[1,1]				# [yr cm^2]
  StErrZeta = header[1,2]			# [yr cm^2]
  REZeta = StErrZeta / Zeta		# [%]
  SqSize = NULL 			# [cm^2]
  
  RhoD = abs (header[2,1])		# [tracks/cm^2]
  StErrRhoD = header[2,2]
  RERhoD = StErrRhoD/RhoD		# [%]
  
  # EffectiveUContent = header[1,3]	# [ppm]
  
  Zeta0 = Zeta
  REZeta0 = REZeta
  RhoD0 = RhoD
  RERhoD0 = RERhoD
  
  K = (Zeta0 * RhoD0) / (Zeta * RhoD)
  # SqNum is known as Dpar in Vers radialplotter
  inputFileData = list(file=filename, details=details, nS = data[,1] , nI = data[,2], SqNum = -1, RhoD = RhoD, Zeta  = Zeta, REZeta = REZeta , RERhoD = RERhoD , SqSize = SqSize , StErrZeta = StErrZeta, K=K, StErrRhoD= StErrRhoD)
  
  return(inputFileData)
}

# # # # # # # # # # # # # # # # # # # # # # # # # 


readFTZFile <- function(dir, filename, ending=".FTZ", seperatorHeader=",", seperatorData=""){
  
  
  tmp <- file(paste(dir,filename, ending ,sep=""), "rt")   ;  details = readLines(tmp, 1)  ;  close(tmp)
  header <- read.table( paste(dir,filename, ending, sep=""), skip=1, sep= seperatorHeader, nrows=2)
  data <- read.table( paste(dir,filename, ending, sep=""), skip=3, sep= seperatorData)
  
  if(header[1,1]==-1){
    # Use Z method
    RERhoD = 0
    print("WARNING:: Z method!!!")
    inputFileData = NULL
    
    SqSize = header[1,3] 
  } else {
    print("Using Zeta method")
    # Use Zeta method
    Zeta = header[2,1]				# [yr cm^2]
    StErrZeta = header[2,2]			# [yr cm^2]
    REZeta = StErrZeta / Zeta		# [%]
    SqSize = header[2,3] 			# [cm^2]
    
    RhoD = abs (header[1,1])		# [tracks/cm^2]
    StErrRhoD = header[1,2]
    RERhoD = StErrRhoD/100		# [%]
    
    EffectiveUContent = header[1,3]	# [ppm]
    
    Zeta0 = Zeta
    REZeta0 = REZeta
    RhoD0 = RhoD
    RERhoD0 = RERhoD
    
    K = (Zeta0 * RhoD0) / (Zeta * RhoD)
    inputFileData = list(file=filename, details=details, nS = data[,1] , nI = data[,2], SqNum = data[,3], RhoD = RhoD, Zeta  = Zeta, REZeta = REZeta , RERhoD = RERhoD , SqSize = SqSize , StErrZeta = StErrZeta, K=K)
  }
  
  return(inputFileData)
}



readAlpsFile <- function( i=2 ){
  
  dataDir = "/Users/mnaylor/Documents/MyPapers/PaperSubmitted/SinclairAndNaylor_2010_bias/AlpineRawData/"
  dataFiles=c("99-11","C7-2(1)","E2-1A(1)","N13-1(1)","N6-2(1)","99-12","C9-2(1)","E2-1C(1)","N13-3(1)","99-22","E1-2(1)","E2-2A(1)","N5-2(1)","C10-2(1)","E18-2(1)","N12-1B(1)","N5-4(1)")
  
  file = dataFiles[i]
  inputData = readBinimfitOutputAsCSV(dataDir, file)
  return(inputData)
}

readBinimfitOutputAsCSV <- function(dataDir, dataFilename){
  filename = paste( dataDir, dataFilename, '.csv', sep="")
  
  paras = scan(filename, skip=10, nlines=1, what=list(x="",y=0, z=""), sep=",")
  ZetaType = paras$x[1]    
  Zeta = paras$y[1]  
  
  paras = scan(filename, skip=11, nlines=1, what=list(x="",y=0, z="", z="", z="", z="", z="", z="", z="" ), sep=",")    
  StErrZeta = paras$y[1]  
  REZeta = StErrZeta / Zeta		# [%]
  
  paras = scan(filename, skip=12, nlines=1, what=list(x="",y=0, z="", z="", z="", z="", z="", z="", z="" ), sep=",")    
  nGrains     = paras$y[1]
  
  paras = scan(filename, skip=16, nlines=1, what=list(x="",y=0, z="", z="", z="", z="", z="", z="", z="" ), sep=",")    
  Nd     = paras$y[1]
  
  paras = scan(filename, skip=17, nlines=1, what=list(x="",y=0, z="", z="", z="", z="", z="", z="", z="" ), sep=",")    
  RhoD    = paras$y[1]
  
  RERhoD = 0
  
  K=1
  
  data = read.csv(file=filename, skip=22, header=TRUE, sep=",", nrows=nGrains)
  
  # Vector of radioactive contect via a proxy nI
  nI = data$Ni
  nS = data$Ns
  grainAges = data$age
  inputFileData = list(nI = nI, nS = nS, grainAges = grainAges, RhoD = RhoD, RERhoD = RERhoD, Zeta  = Zeta, REZeta=REZeta ,StErrZeta = StErrZeta, Nd=Nd , K=K )
  return(inputFileData)
}
