"crosshyb_func" <-
function(input, seqList, mysub, numPermut, probes, satValue, maxProbes, delta=10){

  # Extract sequences
  seqs <- sapply(seqList, function(z){as(z, "character")})

  # To avoid Inf when computing log, replace 0 by 1
  input$rawValue[input$rawValue == 0] <- 1

  # Permutation test to compute p.values
  # The test statistic is the sum of putative corrupted raw values
  perm    <- list()
  refSet  <- input$rawValue[input$maControls %in% probes & !duplicated(input$Name)] # Remove duplicates
  mapping <- c(1:length(input$Name))[input$maControls %in% probes]

  # Sort by raw value,
  ordByValue <- order(input$rawValue, decreasing=TRUE)

  # Keep selected probes (default = c("probes", "spikes")
  ordByValue <- ordByValue[input$maControls[ordByValue] %in% probes]

  # Find the number of genes to consider using the maximum between:
  # - the number of genes whose log2(raw value) is > log2(mean) + 3sd
  # - the number of genes with raw value > satValue (default: 65535)
  sdThr        <- 2^(mean(log2(input$rawValue[ordByValue]), trim=0.1) + 3*sd(log2(input$rawValue[ordByValue])))
  parentProbes <- ordByValue[input$rawValue[ordByValue] >= min(sdThr, satValue)]

  # If two or more 'parents' have the same sequence,
  # keep only that with the highest raw values
  parentProbes <- (parentProbes)[!duplicated(seqs[ match(parentProbes, mapping) ])]

  # Keep a maximun of 'maxProbes' parentProbes to analyze
  if (length(parentProbes) > maxProbes){ parentProbes <- parentProbes[1:maxProbes] }

  # Convert probe number to probe name
  parentNames  <- input$Name[parentProbes]

  print(paste("Putative corruptors to analyze:", length(parentProbes)))

  ptest    <- vector()
  children <- list()

  for (i in 1:length(parentProbes)){
    print(i)
    # Find similar probes (children) based on TagName
    # childrenProbes <- unique(c(grep(substr(parentNames[i], 1,3), input$Name),
    #                           grep(substr(parentNames[i], 2,4), input$Name)))

    # Find similar probes (children) based on Delta Tm < delta
    parentSeq <- seqList[[ match(parentProbes[i], mapping) ]]
    parentTm  <- salign(parentSeq, parentSeq, mysub, alignment="local", delta=-20, gapext=-10, scoring="score")
    allTm     <- salign(parentSeq, seqList, mysub, alignment="local", delta=-20, gapext=-10, scoring="score")
    deltaTm   <- parentTm - allTm
    childrenProbes <- mapping[deltaTm <= delta]

    # Remove 'corruptor' from the 'corrupted' list
    #childrenProbes <- childrenProbes[! input$Name[childrenProbes] %in% parentNames[i]]

    # Remove 'corrupted' whose sequence == 'corruptor' sequence (including the 'corruptor' itself)
    childrenProbes <- childrenProbes[! seqs[match(childrenProbes, mapping)] %in% as(parentSeq, "character")]

    # Keep 'children' with raw values < parent raw value
    childrenProbes <- childrenProbes[ input$rawValue[childrenProbes] < input$rawValue[parentProbes[i]] ]

    # If two or more 'children' have the same sequence,
    # keep only that with the highest raw values
    ord <- order(input$rawValue[childrenProbes], decreasing=TRUE)
    childrenProbes <- (childrenProbes[ord])[!duplicated(seqs[ match(childrenProbes[ord], mapping) ])]
  
    # Update list
    children[[i]] <- childrenProbes

    # Permutation test
    if (length(childrenProbes) > 0){
      childrenValues <- input$rawValue[childrenProbes]
      len <- length(childrenValues)
      if (is.null(perm[[paste("size",len,sep="")]])){
        for (j in 1:numPermut){
          perm[[paste("size",len,sep="")]][j] <- sum(sample(refSet, len))
        }
      }
      ptest[i]   <-  sum(perm[[paste("size",len,sep="")]] > sum(childrenValues))/numPermut
    }else{
      ptest[i]   <- NA
    }
  }

  # Replace 0 values
  ptest[ptest == 0] <- 1/(numPermut * 10)

  # Multiple test correction
  ptest.adj <- p.adjust(ptest, method="BH")

  result <- list(data=data.frame( parentProbes,
                                  parentNames,
                                  ptest.adj,
                                  row.names=parentProbes
                                 ),
                children=children)
  return(result)
}

