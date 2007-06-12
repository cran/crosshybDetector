##Copyright 2001, W. Wolski, all rights reserved.
##
# Class Submatrix
##
library("methods") 

setClass("Submatrix",
         representation(copyright="character",info="character",head="character",alphabet="character")
         ,contains="matrix",prototype(copyright="GNU GENERAL PUBLIC LICENSE Version 2, June 1991")
         )

#promptClass("Submatrix")

if (!isGeneric("subFromEmboss"))
    setGeneric("subFromEmboss",
               function(object,path,...)
               standardGeneric("subFromEmboss"))



setMethod("subFromEmboss",signature(object="Submatrix",path="character"),
          function(object,path,...)
          {
            con<-file(path,"r")
            res <- readLines(con=con,n=-1)
            close(con)
            path<-unlist(strsplit(path,"/"))
            object@info<-path[length(path)]
            head <- res[grep("#",res)]
            head<-paste(head,collapse="\n")
            head<-paste(head,"\n",sep="")
            object@head <- head
            res <- res[-grep("#",res)]
            alphabet <- unlist(strsplit(res[1]," +"))
            object@alphabet <- alphabet[2:length(alphabet)]
            res2<-NULL
            for(x in 2:length(res))
              {
                if(nchar(res[x])!=0)
                  {
                    tmp<-unlist(strsplit(res[x]," +"))
                    res2<-rbind(res2,as.numeric(tmp[2:length(tmp)]))
                  }
              }
            object@.Data <- res2
            colnames(object)<-object@alphabet
            rownames(object)<-object@alphabet
            object
          })

#promptMethods("subFromEmboss")
setMethod("show","Submatrix",function(object)
          {
            cat("info : ",object@info,"\n")
            cat("copyright : ",object@copyright,"\n")
            cat("head :\n" ,object@head)
            cat("alphabet : ",paste(object@alphabet,collapse=" "),"\n")
            print(as(object,"matrix"))
          })

Submatrix<-function()
{
  new("Submatrix")
}

#
#class AAAlphabet
#
setClass("AAAlphabet",representation(info="character")
         ,contains="character",
         prototype(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X","*")
                   ,info="AminoAcid"
                   )
         )

#promptClass("AAAlphabet")

setMethod("levels",signature(x="AAAlphabet"),
          function(x)
          {
            return(levels(as.factor(unlist(strsplit(x,"")))))
          }
          )
setMethod("show",signature(object="AAAlphabet"),
          function(object)
          {
            cat("info : ",object@info,"\n")
            cat("Alphabet :\n",paste(as(object,"character"),collapse=" "),"\n")
          }
          )

##
# AASequence
##

setClass("AASequence"
         ,representation(
                         info="character"
                         )
         ,contains="character"
         )


#promptClass("AASequence")
AASequence <- function(info,sequence)
  {
    if(!missing(info) && !missing(sequence))
      new("AASequence",sequence,info)
    else if(!missing(sequence))
      new("AASequence",sequence)
    else
      new("AASequence")
  }

setMethod("initialize"
          ,signature(.Object="AASequence")
          ,function(.Object,sequence,info,alphabet=new("AAAlphabet"))
          {
            if(!missing(sequence))
              {
                sequence<-toupper(sequence)
                                        #first check if the string is constructed out of alphabet letters.
                seqlevel <- as.character(levels(as.factor(unlist(strsplit(sequence , "" )))))
                alphlevel <- levels(alphabet)
                if(sum(is.element(seqlevel,alphlevel)) != length(seqlevel))
                  {
                    print(sequence)
                    stop("This chars are not in the alphabet: ",paste(seqlevel[is.element(seqlevel,alphlevel)==FALSE]
                                                                      ,collapse=" "
                                                                      ,"!!!\n"))
                  }
                .Object@.Data <- sequence
              }
            if(!missing(info))
              .Object@info<-info
            .Object
          }
          )


setMethod("show",signature(object="AASequence")
          ,function(object)
          {
            cat("info : ", object@info,"\n")
            if(nchar(object)>23)
              cat("sequence:\n",substr(object,1,10),"...",substr(object,nchar(object)-10,nchar(object)),"\n",sep="")
            else
              cat("sequence:\n",as(object,"character"),"\n")
          }
          )



setMethod("frequency",signature(x="AASequence"),
          function(x,alphabet=new("AAAlphabet"),...)
          {
            res<-table(strsplit(x,""))
            res <- res[alphabet]
            res[is.na(res)]<-0
            res<-as.numeric(res)
            names(res)<-alphabet
            return(res)
          }
          )




##
## computes the score of the sequence with itself.
##

if (!isGeneric("selfalign"))
    setGeneric("selfalign",
               function(object,sub,...)
               standardGeneric("selfalign"))

setMethod("selfalign"
          ,signature(object="AASequence",sub="Submatrix")
          ,function(object,sub)
          {
            sdiag<-sub[diag(dim(sub)[1])==1]
            names(sdiag)<-colnames(sub)
            fre<-frequency(object)
            res<-NULL
            res<-0
            for(x in names(fre))
              {
                res<-res+ sdiag[x]*fre[x]
              }
            return(res)
          }
          )

#promptMethods("selfalign")


setClass("AAAlignment"
         ,representation(
                         info1="character"
                         ,info2="character"
                         ,selfs1="numeric"
                         ,selfs2="numeric"
                         ,score="numeric"
                         ,identity="numeric"
                         ,alignsimilarity="numeric"
                         ,lch1="numeric"
                         ,lch2="numeric"
                         ,alig1="character"
                         ,alig2="character"
                         ,beautify="character"
                         )
         )

#promptClass("AAAlignment")
          ##
          ##v \item{similarity}{
          ##v This is a count of the number of positions over the length of the alignment where >= 51% of the residues or bases at that position are similar. 
          ##v Any two residues or bases are defined as similar when they have positive comparisons (as defined by the comparison matrix being used in the alignment algorithm). 
          ##v  }
          ##v \item{identity}{
          ##v This is a count of the number of positions over the length of the alignment where all of the residues or bases at that position are identical. 
          ##v }

setMethod("show"
          ,signature(object="AAAlignment")
          ,function(object)
          {
            lalign<-nchar(object@alig1)
            cat("selfscore 1: ",object@selfs1,"; seq length 1 :",object@lch1,"\n")
            cat("selfscore 2: ",object@selfs2,"; seq length 2 :",object@lch2,"\n")
            cat("alig lenght: ",lalign,"\n")
            cat("score      : ",object@score,"\n")
            cat("FM         : ",object@score/(sqrt(object@selfs1*object@selfs2)),"\n")
            cat("identity   : ",object@identity,"/",min(object@lch1,object@lch2),"\n")
            cat("similarity : ",object@alignsimilarity,"/",min(object@lch1,object@lch2),"\n")
          }
          )

setMethod("summary"
          ,signature(object="AAAlignment")
          ,function(object)
          {
            identity <- sum(unlist(strsplit(object@alig1,""))==unlist(strsplit(object@alig2,"")))
            nammax <- max(nchar(object@info1),nchar(object@info2))
            lalign<-nchar(object@alig1)
            cat("selfscore 1: ",object@selfs1,"\n")
            cat("selfscore 2: ",object@selfs2,"\n")
            cat("alig lenght: ",lalign,"\n")
            cat("score      : ",object@score,"\n")
            cat("FM(score)  : ",object@score/(sqrt(object@selfs1*object@selfs2)),"\n")
            cat("identity   : ",object@identity,"/",min(object@lch1,object@lch2),"\n")
            cat("similarity : ",object@alignsimilarity,"/",min(object@lch1,object@lch2),"\n")
            tmp<-c(seq(1,lalign,60),lalign)
            for(x in 1:(length(tmp)-1))
              {
                s1 <- substr(object@alig1 , tmp[x] , tmp[x+1] )
                s2 <- substr(object@alig2 , tmp[x] , tmp[x+1] )
                beauti <- substr(object@beautify ,tmp[x],tmp[x+1])
                cat( format( object@info1 , width=nammax ) , s1 , "\n" )
                cat( format( " ", width =nammax), beauti, "\n")
                cat( format( object@info2 , width=nammax ) , s2 , "\n" )
                cat("\n")
              }
          }
          )

##
##Class : AASequenceList
##

setClass("AASequenceList"
         ,representation(info="character",names="character"
                     #    ,list="list" #change
                         )
         ,contains="list"
         )


setMethod("show","AASequenceList"
	,function(object)
	{
		cat("info:    ", object@info,"\n")
		cat("length : ", length(as(object,"list")),"\n")
	}
	)	


setReplaceMethod("[[", "AASequenceList"
                 , function(x, i, j,..., value)
                 {
                   if( !extends(class(value),"AASequence") )
                     {
                       stop(paste("This is an AASequenceList!"
                                  ,"so dont try to assing a object of class:\",the object is class"
                                  ,class(value)
                                  ,"\n"
                                  ,sep=" ")
                            )
                     }
                   tt<-as(x,"list")
                   tt[[i]]<-value
                   names(tt)[i]<-value@info
                   as(x,"list")<-tt
                   x
                 })

          

setMethod("[",
          "AASequenceList",
          def = function(x, i, j, ..., drop = F)
          {
            y <- as(x,"list")
            names(y)<-names(x)
            as(x,"list") <- y[i]
            return(x)
          }
          )

setReplaceMethod("[","AASequenceList"
                 ,function(x,i,j,...,value)
                 {
                   if( !extends(class(value),"AASequenceList") )
                     {
                       stop(paste("This is an AASequenceList!"
                                  ,"so dont try to assing a object of class:\n"
                                  ,class(value)
                                  ,"\n"
                                  ,"Only Objects of class AASequenceList can be assigned.\n"
                                  ,sep=" ")
                            )
                     }
                   y<-as(x,"list")
                   y[i] <- value
                   as(x,"list") <- y
                   x
                 }
                 )

##
##Alignment Methods
##

if (!isGeneric("salign"))
    setGeneric("salign",
               function(obj1,obj2,...)
               standardGeneric("salign")
               )

setMethod("salign"
          ,signature(obj1="AASequenceList",obj2="NULL")
          ,function(obj1
                    ,obj2
                    ,sub
                    ,delta= -4
                    ,gapext = delta
                    ,alignment = "global"
                    ,scoring = "identity"
                    ,diag = FALSE
                    )
          {
            res<-listdist(obj1
                          ,testalign
                          ,diag=diag
                          ,sub
                          ,delta=delta
                          ,gapext=gapext
                          ,alignment=alignment
                          ,scoring=scoring
                          )
            if(scoring=="identity"||scoring=="similarity"||scoring=="scoreN")
              {
                res[1:length(res)] <- (1-as.numeric(res))
              }
            else if(scoring=="score")
              {
                tmp <- as.numeric(res)
                tmp <- (tmp - mean(tmp))/sqrt(var(tmp)) 
                res[1:length(res)] <- pnorm(tmp,mean=0,sd=1,lower.tail=FALSE)
              }
            else if(scoring=="pozitive")
              {
                res[1:length(res)] <- pnorm(res,mean=0,sd=1,lower.tail=FALSE)
              }
            return(res)
          }
          )

if (!isGeneric("listdist"))
    setGeneric("listdist",
               function(object,...)
               standardGeneric("listdist"))



setMethod("listdist"
          ,signature(object="list")
          ,function(object,FUN,diag=FALSE,...)
          {
            lo<-length(object)
            if(length(object)==1)
              return(dist(1))
            res<-numeric(lo*(lo-1)/2)
            aa <- 1
            for(rr in 1:lo)
              {
                tt<-(rr+1)
                if(tt <= lo)
                  {
                    SL <- object[tt:lo] #sublist
                    ee <- (aa-1) + length(SL)
                    tmp <- object[[rr]]
                    res[aa:ee] <- unlist(lapply(SL,FUN,tmp,...))
                  }
                aa <- ee + 1
              }
            ans <- res
            attributes(ans) <- NULL
	    attr(ans,"Labels") <- names(object)
            attr(ans,"Size") <- length(object)
            attr(ans, "call") <- match.call()
            class(ans) <- "dist"
            attr(ans,"Diag") <- diag
            attr(ans,"Upper") <- TRUE
            return(ans)
          }
          )


setMethod("salign"
          ,signature(obj1="AASequence",obj2="AASequenceList")
          ,function(obj1,obj2,sub,delta=-4,gapext = delta, alignment="global",scoring="pozitive")
          {
            res<-salign(obj2,obj1,sub, delta=delta , gapext=gapext , alignment=alignment,scoring = scoring)
          }
          )

setMethod("salign"
          ,signature(obj1="AASequenceList",obj2="AASequence")
          ,function(obj1,obj2,sub,delta=-4,gapext = delta, alignment="global",scoring="pozitive")
          {
            res<-lapply(obj1,testalign,obj2,sub,delta=delta,gapext=gapext,alignment=alignment,scoring=scoring)
            res <- unlist(res)
          }
          )


##t Pairwise sequence Aligment of Amino Acid Sequence
##- Pairwise sequence Aligment of Amino Acid Sequence. Function can compute
##- global, local or overlap alignment of two amino acid sequences.
##+ ret : what to return? AAAlignment = object of class AAAlignement, identity/alignment length, similarity, score.
##+ delta : gap opening penalty
##+ gepext : gap extension penalty
##+ type :  type of alignemnt. (e.g. global,local,overlap)
##+ sub : similarity matrix BLOSUM; PAM similarity matrix.
##+ obj1 : object of class AASequence
##+ obj2 : object of class AASequence
##e 
##e 

setMethod("salign",signature(
                            obj1="AASequence"
                            ,obj2="AASequence"
                            )
          ,function(
                    obj1
                    ,obj2
                    ,sub
                    ,delta=-4
                    ,gapext = delta
                    ,alignment="global"
                    ,scoring="AAAlignment"
                    )
          {
            mret<-c("AAAlignment","identity","similarity","score","scoreN")
            if(!(scoring %in% mret))
              {
                stop("scoring argument can be either ", paste(mret,collapse=" "),"\n")
              }
            
            ####real computing.
            res<-.Call("alignSEXP"
                       ,obj1
                       ,obj2
                       ,sub
                       ,delta
                       ,gapext
                       ,alignment
                       )

            if(scoring=="AAAlignment")
              {
                 if(length(grep("   ",res[["errmsg"]])) == 0)
                   {
                     stop("ERROR:", errmsg, "\n")
                   }
           
                ss1<-paste(res[["alig1"]],"*",sep="")
                ss2<-paste(res[["alig2"]],"*",sep="")
                lalign<-nchar(ss1)
                tmp<-c(seq(1,lalign,40),lalign)
                wmatch<-rep(" ",lalign)
                beauti<-""
                for(x in 1:(length(tmp)-1))
                  {
                    vs1 <- substr(ss1 , tmp[x] , tmp[x+1]-1 )
                    vs2 <- substr(ss2 , tmp[x] , tmp[x+1]-1 )
                    vs1<-unlist(strsplit(vs1,""))
                    vs2<-unlist(strsplit(vs2,""))
                    match <- wmatch[tmp[x]:(tmp[x+1]-1)]
                    ##find similarities.
                    ## what exactly means the star in the blosum matrix.
                    ## you anyway will count only positive values
                    ## * AA are always smaller than 0 so it does what we are looking for.
                    vs1t<-vs1
                    vs2t<-vs2
                    vs1t[vs1=="-"] <- "*"
                    vs2t[vs2=="-"] <- "*"
                    tmpsim <- sub[vs1t,vs2t] # get the values of the diagonal
                    
                    tmpsim <- tmpsim[diag(dim(tmpsim)[1])==1] # get the diagonal
                    names(tmpsim)<-NULL
                    match[tmpsim>0]<-":" #mark similarities
                    match[vs1==vs2]<-"|" #mark identities
                    beauti<-paste(beauti,paste(match,collapse=""),sep="")
                  }
                res<- new("AAAlignment"
                          ,info1 = obj1@info
                          ,info2 = obj2@info
                          ,selfs1 = res[["selfscore1"]]
                          ,selfs2 = res[["selfscore2"]]
                          ,score=res[["score"]]
                          ,identity = res[["identity"]]
                          ,alignsimilarity = res[["alignsimilarity"]]
                          ,alig1 = res[["alig1"]]                      
                          ,alig2 = res[["alig2"]]
                          ,lch1= nchar(obj1)
                          ,lch2= nchar(obj2)
                          ,beautify = beauti
                          )
                return(res)
              }
            else if(scoring=="similarity")
              {
                return(res[["alignsimilarity"]]/min(nchar(obj1),nchar(obj2)))
              }
            else if(scoring=="identity")
              {
                return(res[["identity"]]/min(nchar(obj1),nchar(obj2)))
              }
            else if(scoring=="scoreN")
              {
                sc <- ifelse(res[["score"]]>0,res[["score"]],0)
                return(sc/min(res[["selfscore1"]],res[["selfscore2"]]))
              }
            else if(scoring=="score")
              {
                return(res[["score"]])
              }
          }
          )
#promptMethods("align")

##
##Read write functions for Fasta format
##
if (!isGeneric("testalign"))
    setGeneric("testalign",
               function(obj1,obj2,...)
               standardGeneric("testalign")
               )

setMethod("testalign"
          ,signature(obj1="AASequence",obj2="AASequence")
          ,function(
                    obj1
                    ,obj2
                    ,sub
                    ,delta=-4
                    ,gapext = delta
                    ,alignment="global"
                    ,scoring="score"
                    )
          {
            res<-.Call("alignScoreSEXP"
                       ,obj1
                       ,obj2
                       ,sub
                       ,delta
                       ,gapext
                       ,alignment
                       ,scoring
                       )
            return(res)
          }
          )



if (!isGeneric("readFasta"))
    setGeneric("readFasta",
               function(object,file,...)
               standardGeneric("readFasta"))

infogrep <- function(x)
  {
    return(sub("^>([a-zA-Z0-9]+) .+","\\1",x,perl=TRUE))
  }

seqgrep <- function(x)
  {
     return(gsub("\\*","",x))
  }

setMethod("readFasta"
          ,signature(object="AASequenceList")
          ,function(object
                    ,file
                    ,grepinfo=infogrep
                    ,grepseq=seqgrep)
          {
            con <-file(file,"r" )
            all <- readLines(con,n=-1)
            pos <-  grep(">",all)
            dat <- vector("list",(length(pos)-1))
            nam <- character(length(pos)-1)
            
            if(length(pos)>1)
              {
                for(x in 1:(length(pos)-1))
                  {
                    info <- grepinfo(all[pos[x]]) # get the info
                                        #cat("x = ",x, "  | pos[x] = ",pos[x], " | ", info , "\n" )
                    seq <- paste(all[(pos[x]+1):(pos[x+1]-1)],collapse="")
                    seq <- grepseq(seq)
                    tmp <- AASequence(info,seq)
                    nam[x] <- info
                    dat[[x]] <- tmp
                  }
                  names(dat)<-nam
              }
	      else
            {	
              
            }
            as(object,"list") <- dat
            return(object)
          }
          )

#promptMethods("readFasta")




