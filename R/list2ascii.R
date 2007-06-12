"list2ascii" <-
function(x,file) {

   # Write an R list to an ASCII file.
   # ASCII files can be used by 
   # a spreadsheet or other programs.
   #
   tmp.wid = getOption("width")  # save current width
   options(width=10000)          # increase output width
   sink(file)                    # redirect output to file
   print(x)                      # print the object
   sink()                        # cancel redirection
   options(width=tmp.wid)        # restore linewidth
   return(invisible(NULL))       # return (nothing) from function

}

