# Appendix: Analysis Functions

## The following data analysis methods were included to clarify data-processing, and to provide direction to others in completing analysis. 

# 1.	Load Packages
library(writexl) # for writing dataframe to excel
library(signal) # for filter function

# 2.	Load Data
## create a dataframe of  names for lubricant group and compression
read_In <- function(date, lubricant, compression, samples){
	
	samplesize <- samples
	dfnamesFT <- c( )
	dfnamesSR <- c( )
	df <- c( )
	df2 <- c( )

  	for (i in 1:samplesize){
		dfnamesFT <- paste0(date,"_TEC_", lubricant, "_S", i, "_", compression, "_FT.txt")
    		dfnamesSR <- paste0(date,"_TEC_", lubricant, "_S", i, "_", compression, "_SR.txt")
    		df2 <- cbind(dfnamesFT,dfnamesSR)
    		df <- rbind(df, df2) 
  		}
  
	return(df)
 }

## read-in friction test file
read_FT <- function(name){
	
	df <- read.delim(name, skip = 22, stringsAsFactors = FALSE)
  	
	# convert any data that are not numbers to numbers or NAs
 	df$Time..s <- as.numeric(as.character(df$Time..s))
 	
	## get rid of rows with NAs or missing data
 	df <- na.omit(df)
 
	return(df)
}

## read-in stress-relaxation file 
read_SR <- function(name){
 	
	df <- read.delim(name, skip = 24, stringsAsFactors = FALSE)
 	df$Time..s <- as.numeric(as.character(df$Time..s))
	df <- na.omit(df)

	return(df)
}

# 3. Pre-Processing Data – Setting Filter
## set filter
bf <- butter(n = 1, W= 0.11, type = "low", plane = "z") 


# 4. Compressive Stiffness
## functions to prepare data for compressive stiffness calculation 
StressStrain <- function(df, COMP){
 
## initialize local variables that will be the output of this function
	Stress <- 0
	Strain <- 0
	Time <- 0
	Force <- 0

## surface area based on scaffold radius, 0.005m
  	Area <- 3.14*(0.005*0.005)
  
## remove early/end data where measurement fluctuations are high
 	MaxTime <- max(df$Time..s) - 0.02
	df <- df[df$Time..s > 0.02 & df$Time..s < MaxTime, ]
  	
## setting local variables from dataframe
	time <- df$Time..s
  	displacement <- df$Position..z...mm
  	length0 <- displacement[5]  
  
## “zero” force 
	Fnot <- df[df$Time..s < 0.5, "Fz..N"]
	FNot <- mean(Fnot)
	Force <- filter(bf, df$Fz..N) - FNot

## setting break value based on compressive stiffness input
	Comp <- COMP/100 - 0.001

## generating stress, strain, force and time dataframe that will be output from this function

	for (i in 5:length(Force)){
    		Stress[i] <- Force[i]/Area
    		Strain[i] <- abs((length0 - displacement[i])/length0)
    		Force[i] <- Force[i]
    		Time[i] <- time[i]
   			if (Strain[i] > Comp) {
     			break
    			}
 	 	}
  
	result <- cbind.data.frame(Time, Force, Strain, Stress)
  
	}

## calculating compressive stiffness from stress-strain data frame 

SR_Stiffness_Function <- function(df){
	
	# filtering values between strain = 0.02 and 0.09, to create new df for analysis

	stiffness <- df[df$Strain > 0.02 & df$Strain < 0.09,]

	# completing simple linear regression, lm function, to determine compressive stiffness 
	reg <- lm(Stress~Strain, data = stiffness)  

	# returning linear regression ß1 coefficient 
	stiff <- reg$coefficients[2]
  
	result <- cbind.data.frame(stiff)  
  
	}

#5. Calculation Normal Equilibrium Load After Stress Relaxation

## preparing file for Neq calculation

NormalEquilibrium <- function(df){
  	
	MaxTime <- max(df$Time..s) - 0.02
	df <- df[df$Time..s > 0.02 & df$Time..s < MaxTime, ]
  	time <- df$Time..s
  
  	Fnot <- df[df$Time..s < 1.5, "Fz..N"]
	FNot <- mean(Fnot)
	Force <- filter(bf, df$Fz..N) - FNot
  
  result <- cbind.data.frame(time, Force)
  
}

## calculating normal equilibrium load, Neq

SR_Neq_Function <- function(df){
  	
	#Neq calculated from the last 10 seconds of data
 	NeqTime <- max(df$time) – 10
	NeqCalculation <- df[df$time > NeqTime,]
  
  	Neq <- mean(NeqA$Force)

 	result <- cbind.data.frame(Neq)  
}

# 6. Calculation Torque from Friction Test

## preparing data for friction test analysis

FT_Function <- function (df) {
  
	Torque <- filter(bf, df$Tz..N.mm) 
	Force <- filter(bf, df$Fz..N) 
	Time <- df$Time..s
	Degree <- df$Position..theta...deg
  
	df <- cbind.data.frame(Time, Degree, Force, Torque)
  
}

## calculating torque and normal force during rotation
TorqueFunction <- function (df) {
 
 	time <- max(df$Time)
	halftime <- time/2
  
	# separating CW and CCW rotations
	Positive <- df[df$Time < halftime,]
	Negative <- df[df$Time > halftime,]
	
	PosT <- Positive[Positive$Degree > 100 & Positive$Degree < 260, "Torque"]
	PosTAvg <- mean(PosT)
 	PosF <- Positive[Positive$Degree > 100 & Positive$Degree < 260, "Force"]
	PosFAvg <- mean(PosF100260)

	NegT <- Negative[Negative$Degree > 100 & Negative$Degree < 260, "Torque"]
	NegTAvg <- mean(NegT)
	NegF <- Negative[Negative$Degree > 100 & Negative$Degree < 260, "Force"]
	NegFAvg <- mean(NegF)
	
  
  results <- cbind.data.frame(PosTAvg, NegTAvg, PosFAvg, NegFAvg)
  
}


# 7. Analysis Function 

## analysis function  to automatically return one row of a dataframe. The values include sample No., compression, compressive stiffness, Neq, positive and negative torque, and force
Analysis <- function(FT, SR, name, compression){

	# format data for stiffness calculation, then calculate, returning S
	S1 <- StressStrain(SR,compression)
	S <- SR_Stiffness_Function(S1)

## format data for Neq calculation, then calculate, returning SR
	SR1 <- NormalEquilibrium(SR)
	SR <- SR_Neq_Function(SR1)

## format data for torque calculation, then calculate, returning FT
	FT1 <- FT_Function(FT)
	FT <- TorqueFunction(FT1)

## bind dataframe values and add sample No./ name 
	return <- cbind.data.frame(name, S, SR, FT)

	}



# 8. Example 
## Name files in this format when completing tests
"10302019_TEC_PBS_S8_50_SR”
"date_TEC_lubricant_sample number_compression_SR stress relaxation file or FT torque file ”

## set date, lubricant, and number of samples
date = "10302019"
lubricant = "PBS"
samples = 3

## set % compression
compression = "10"
COMP = 10

## create dataframe of file names
df <- read_In(date, lubricant, compression, samples)

## read in data for samples 1 to 3
S1_FT <- read_FT(df[[1,1]])
S1_SR <- read_SR(df[[1,2]])
S2_FT <- read_FT(df[[2,1]])
S2_SR <- read_SR(df[[2,2]])
S3_FT <- read_FT(df[[3,1]])
S3_SR <- read_SR(df[[3,2]])

## complete analysis on files
R1 <- Analysis(S1_FT, S1_SR, paste0("S1_",compression, "_", lubricant), COMP)
R2 <- Analysis(S2_FT, S2_SR, paste0("S2_",compression, "_", lubricant), COMP)
R3 <- Analysis(S3_FT, S3_SR, paste0("S3_",compression, "_", lubricant), COMP)

S <- rbind.data.frame(R1, R2, R3)

## export results as an excel file
write_xlsx(S, paste0("Summary_",lubricant,"_", compression, "_Results.xlsx"))
