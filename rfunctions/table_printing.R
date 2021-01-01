


print_out_table <- function(x1, x2, x3, x4, x5, x6,....) {
	
	ndec <- function(x, k) {
		trimws(format(round(x, k), nsmall=k))
	    }

	cat("Mortality  &", ndec(x1[1], 2), "&",  ndec(x4[1], 2), " \\\\ 
  & [", ndec(x1[2], 1), ",", ndec(x1[3], 1), "] & [",
    ndec(x4[2], 1), ",", ndec(x4[3], 1),  "] \\\\ \n") 
    
    cat("Complication  &", ndec(x2[1], 1), "&",  ndec(x5[1], 1), " \\\\ 
  & [", ndec(x2[2], 1), ",", ndec(x2[3], 1), "] & [",
    ndec(x5[2], 1), ",", ndec(x5[3], 1),  "] \\\\ \n") 
    
    cat("Length of Stay &", ndec(x3[1], 1), "&",  ndec(x6[1], 1), " \\\\ 
  & [", ndec(x3[2], 1), ",", ndec(x3[3], 1), "] & [",
    ndec(x6[2], 1), ",", ndec(x6[3], 1),  "] \\\\ \n") 
	
}



print_ittsens_table <- function(x1, x2, x3, x4,....) {
	
	ndec <- function(x, k) {
		trimws(format(round(x, k), nsmall=k))
	    }

	cat("Complication  &", ndec(x1, 2), "&",  ndec(x2, 2), " \\\\  \n") 
        
    cat("Length of Stay &", ndec(x3, 2), "&",  ndec(x4, 2), " \\\\  \n") 
	
}



print_sens_table <- function(x1, x2, x3, x4, x5, x6,x7,x8,....) {
	
	ndec <- function(x, k) {
		trimws(format(round(x, k), nsmall=k))
	    }

	cat("Complication  &", ndec(x1, 2), "&",  ndec(x5, 2), " \\\\  \n") 
    
    cat("Complication + Effect Mod. &", ndec(x2, 2), "&",  ndec(x6, 2), " \\\\ \n") 
    
    cat("Length of Stay &", ndec(x3, 2), "&",  ndec(x7, 2), " \\\\  \n") 
    
    cat("Length of Stay + Effect Mod &", ndec(x4, 2), "&",  ndec(x8, 2), " \\\\  \n") 
	
}


print_sens_em_table <- function(x1, x2, x3, x4, x5, x6,....) {
	
	ndec <- function(x, k) {
		trimws(format(round(x, k), nsmall=k))
	    }
    
    cat("Length of Stay + Conventional SE &", ndec(x1, 2), "&",  ndec(x7, 2), " \\\\  \n") 
    
    cat("Length of Stay + Linear Model Se &", ndec(x4, 2), "&",  ndec(x8, 2), " \\\\  \n") 
    
    cat("Length of Stay + Pairs of Pairs &", ndec(x4, 2), "&",  ndec(x8, 2), " \\\\  \n") 
	
}




    
    
    
    
    


     
