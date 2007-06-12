"crosshyb" <-
function(raw, probeSeq, plate = 1, probeNameID="ProbeName", numPermut = 10000, 
                     probes = c("probes", "spike"), satValue = 65535,
                     maxProbes = 50, delta=10){

  # Create Tm matrix
  # Using this matrix the alignment score is the Tm
  # computed using the Basic Tm calculator:
  # Tm = 4°C  x  (number of G’s and C’s) + 2°C  x  (number of A’s and T’s) 
  data(EBLOSUM50)
  mysub <- EBLOSUM50
  mysub@.Data[]      <- 0
  mysub@.Data[1,1]   <- 2
  mysub@.Data[5,5]   <- 4
  mysub@.Data[8,8]   <- 4
  mysub@.Data[17,17] <- 2

  # Create list of 'AASequence' objects
  print ("Creating fasta file...")
  tmpFile <- tempfile()
  cat(paste(paste(">", (maInfo(maGnames(raw))[[probeNameID]])[maControls(raw) %in% probes], sep=""), probeSeq[maControls(raw) %in% probes], sep="\n"), ">",  sep="\n", file=tmpFile)
  seqList <- new("AASequenceList",info="my sequence list")
  seqList <- readFasta(seqList, file=tmpFile, grepinfo=infogrep, grepseq=seqgrep)
  unlink(tmpFile)
  print("Fasta file created!")

  # Analysis of Red Channel
  print("Analyzing RED channel")
  input <- data.frame(rawValue=as.integer(maRf(raw[,plate])),
                      maControls=maControls(raw),
                      Name=I(maInfo(maGnames(raw))[[probeNameID]]))
  class(input$Name) <- "character"
  resultR <- crosshyb_func(input, seqList, mysub, numPermut, probes, satValue, maxProbes, delta)

  # Analysis of Green channel
  print("Analyzing GREEN channel")
  input <- data.frame(rawValue=as.integer(maGf(raw[,plate])),
                      maControls=maControls(raw),
                      Name=I(maInfo(maGnames(raw))[[probeNameID]]))
  class(input$Name) <- "character"
  resultG <- crosshyb_func(input, seqList, mysub, numPermut, probes, satValue, maxProbes, delta)

  # Return output as list
  result <- list(dataR=resultR$data,
                 dataG=resultG$data,
                 childrenR=resultR$children,
                 childrenG=resultG$children)
  return(result)
}

