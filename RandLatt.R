	bs <- function (S_0,k,T,sig,rf,phi) {
		d1 = (log(S_0/k) + (rf + 0.5*sig^2)*T)/( sig*sqrt(T) )
		d2= d1 - sig*sqrt(T)
		optval = S_0*phi*pnorm(phi*d1) - phi*k*exp(-rf*T)*pnorm(d2*phi)
		optval
	}
	
	crBuckets <- function(S,numB) {
		
		spacing <- (max(S) - min(S)) / (mBuckets - 1)
		max(S) - spacing*(numB-1)
	}
	
	St <- function(S,e,t){
		S[,t-1] * exp( (mu - 0.5*sigma^2) * h + sigma*e[,t-1]*sqrt(h) )	 			 	
	} 
	
	K = c(90,95,100,105,110) 
	
	
	S_0 <- 100				##### initial stock value
	
	mu <- 0.10				##### average value of S
	sigma <- 0.20 			##### volatility of S
	
	T <- 1					##### Time until maturity (years) 
	h <- 1/20				##### time step size
	dTime <- T / h			##### number of steps through time
	
	nRandWalks <- 100000			##### number of random walk simulations
	mBuckets <- 300			##### number of stock value buckets

		
	e <- matrix( rnorm( nRandWalks*dTime), nRandWalks, dTime)
	S <- matrix(0, nRandWalks, dTime + 1)
	
	S[,1] <- S_0
		
#######################################################################################################################	

	for (t in 2:(dTime + 1)){
		S[,t] <- St(S,e,t)
		
	 	if ( t == (dTime + 1) ){
	 		spacing <- ( max(S) - min(S) ) / (mBuckets - 1)
	 		Buckets <- crBuckets(S,1:mBuckets)	
	 	}
	}
	
	counts <- matrix(0,mBuckets,dTime+1) 						##### of out-traversals from Bucket(i,t)
	

	Payoff <- matrix(0,mBuckets,dTime+1)
	St <- matrix(0,mBuckets,dTime+1)
	
	Allprobs <- list()
	Allcounter <- list()
	
	for (t in (dTime+1):2 ){ #:(dTime+1)
		
		counter <- matrix(0,mBuckets,mBuckets)
		probs <- matrix(0,mBuckets,mBuckets)
		
		for (m1 in 1:mBuckets){
			
			indices <- which( abs( Buckets[m1] - S[,t-1]) < spacing/2, arr.ind = TRUE)  ##### find the indices of S[,t-1] that are closest to Bucket[m1] 
			for (m2 in 1:mBuckets){
				counter[m1,m2] <- length( which( abs( Buckets[m2] - S[indices,t]) < spacing/2, arr.ind = TRUE))  	##### number of traversals made from Bucket[m1](t-1) to Bucket[m2](t)
				if ( counter[m1,m2] != 0){
					St[m2,t] <- Buckets[m2]
				}
			}
			
			counts[m1,t-1] <- sum(counter[m1,])	  	# of out-traversals from Bucket[m1] (dTime:1)
			probs[m1,] <- counter[m1,] / counts[m1,t-1]
			NaNs <- which(probs[m1,] == 'NaN', arr.ind = TRUE)
			probs[m1,] <- replace(probs[m1,],NaNs,0)	
			
		}
		
		Allprobs[[t-1]] <- probs
	} 
	
	rootBucket <- which(abs(Buckets - S_0) < spacing/2, arr.ind = TRUE)
	St[rootBucket,1] <- Buckets[rootBucket]

	output <- data.frame()
	
	for (i in seq_along(K)){
		k <- K[i]
		
		for (t in (dTime + 1):2){
			
			probs <- Allprobs[[t-1]]
			
			for (m1 in 1:mBuckets){
				if ( t == dTime+1){
					Payoff[,dTime+1] <-pmax(0,St[,dTime+1] - k)			
					if( counts[m1,t-1] != 0){
						Payoff[m1,dTime] <- pmax( St[m1,dTime] - k, exp(-mu*h)*probs[m1,]%*%Payoff[,dTime+1])
					}	
				}
				else if (t > 2 & t < (dTime+1) ){ 
			
					Payoff[m1,t-1] <- pmax( St[m1,t-1] - k, exp(-mu*h)*probs[m1,]%*%Payoff[,t])

				}
				else{
					Payoff[m1,t-1] <- pmax( St[m1,t-1] - k, exp(-mu*h)*probs[m1,]%*%Payoff[,t] )
			
				}
				
			}
		
		}
		
		simCall <- Payoff[rootBucket,1]
		bsCall <- bs(S_0,k,T,sigma,mu,1)
		ratio <- simCall/bsCall
		
		single <- data.frame("Strike Price" = k, "Black-Scholes" = bsCall ,"RandLatt" = simCall,"Price Ratios" = ratio)
		output <- rbind(output,single)
		
	}
	print(output)
	
	
	
		
