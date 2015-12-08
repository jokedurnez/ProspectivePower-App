################################
# NEUROPOWER: CLUSTER FUNCTION #
################################

# Function: cluster
# -----------------
# This function extracts all peaks from a 3D map of
# statistical values (Z or T) above a certain excursion threshold u
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
# 	- zmap: 3D statistical image
#		- u: excursion threshold

cluster <- function(zmap,excursion){

	# enlarge the map first for easier calculation at borders
	#set outside border at value -100
	options(warn=-1)
	dimZ <- dim(zmap)
	Zstar <- array(0,dim=dimZ+c(2,2,2))
	Zstar[2:(dimZ[1]+1),2:(dimZ[2]+1),2:(dimZ[3]+1)] <- zmap
	Zstar[is.na(Zstar)] <- -100

	# empty container of peaks and coordinates
	peaks <- xco <- yco <- zco <- c()

	# loop over all voxels to see if surrounding voxels are higher
	for(m in 2:(dim(Zstar)[1]-1)){
		# in this loop a progression is calculated
		# to show a progression bar in the webtool
		incProgress(1/m,detail=paste("Percentage finished:",100*round(m/(dim(Zstar)[1]-1),digits=2),"%"))
		for(n in 2:(dim(Zstar)[2]-1)){
			for(o in 2:(dim(Zstar)[3]-1)){
				if(Zstar[m,n,o]<excursion){next} # if the value doesn't exceed
																 # the excursion threshold: next voxel
				surroundings <- c(Zstar[m-1,n+1,o],
				Zstar[m,n+1,o],
				Zstar[m+1,n+1,o],
				Zstar[m-1,n,o],
				Zstar[m+1,n,o],
				Zstar[m-1,n-1,o],
				Zstar[m,n-1,o],
				Zstar[m+1,n-1,o],

				Zstar[m-1,n+1,o-1],
				Zstar[m,n+1,o-1],
				Zstar[m+1,n+1,o-1],
				Zstar[m-1,n,o-1],
				Zstar[m+1,n,o-1],
				Zstar[m-1,n-1,o-1],
				Zstar[m,n-1,o-1],
				Zstar[m,n,o-1],
				Zstar[m+1,n-1,o-1],

				Zstar[m-1,n+1,o+1],
				Zstar[m,n+1,o+1],
				Zstar[m+1,n+1,o+1],
				Zstar[m-1,n,o+1],
				Zstar[m+1,n,o+1],
				Zstar[m-1,n-1,o+1],
				Zstar[m,n,o+1],
				Zstar[m,n-1,o+1],
				Zstar[m+1,n-1,o+1])

				# if the voxel value is higher than all
				# surrounding voxel values: add to peak table
				if(Zstar[m,n,o]>max(surroundings)){
					peaks <- c(peaks,Zstar[m,n,o])
					xco <- c(xco,m-1)
					yco <- c(yco,n-1)
					zco <- c(zco,o-1)
				}
			}
		}
	}
	results <- data.frame(peaks,xco,yco,zco)
	return(results)
}
