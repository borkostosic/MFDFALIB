dirname <- "D:\\laptop_07_02_2010\\D\\papers\\MDFA\\MFDFALIB\\R"
setwd(dirname)

fname <- "..\\DATA\\ser16.txt"		# Binomial multifractal model 2^16=65536 dta poits
seq <- read.table(fname)		# load data

MAX_BOX<-200		# maximum 200 points on logarithmic scale...
MAXQ<-201		# max q resolution -10,...,10 dq=0.1

data<-as.numeric(unlist(seq))
total <- length(data)	# data length
qmin <- -10.000001	# small offset to avoid 0
qmax <- 10.0
dq <- 0.1
qseq <- seq(qmin, qmax, by = dq)	# sequence of q's
nq<-length(qseq)			# length of q sequence
H <- numeric(length(qseq))		# generalized Hurst exponet
tau <- numeric(length(qseq))		# Renyi exponet
f <- numeric(length(qseq))		# f(alpha)
alpha <- numeric(length(qseq))		#alpha
dmse<-as.numeric(matrix(0, nrow=MAXQ, ncol=MAX_BOX))	# fluctuations
nfit<-1		# polynomial degree 1 (linear), 2 or 3
sw<-0		# sw=1 for sliding segments, takes longer, use with care... 

rs<-integer(MAX_BOX)			# segment size list
# prepare multiplicative scale (equidistant on log-log plot)
   nrs<-0
   minseg<-4		# minimum segment size
   maxseg<-total/4	# maximum segment size
   boxratio<-2^(1/8)	# segment ratio (multiplicative scale)
#...OR prepare your own segment scale
   #rs<-10:100 	
   #nrs<-length(rs)

dllname <- paste(dirname,"\\mfdfa64.dll",sep="")	# dll name
dyn.load(dllname)					# load the dynamic library
start.time <- Sys.time()	# call mfdfa_calc...
ans <- .C("mfdfa_calc", as.numeric(data), as.integer(total),
	as.numeric(H), as.numeric(tau), 	
	as.numeric(f), as.numeric(alpha),
	as.numeric(dmse), as.integer(rs), as.integer(nrs),  
	as.numeric(qmin), as.numeric(qmax), as.numeric(dq),
	as.integer(minseg), as.integer(maxseg), as.numeric(boxratio),
	as.integer(nfit), as.integer(sw)
	)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

dyn.unload(dllname)					# load the dll

H<-ans[[3]]
tau<-ans[[4]]
f<-ans[[5]]
alpha<-ans[[6]]
dmse<-ans[[7]]
rs<-ans[[8]]
#nr<-sum(rs != 0)
nr<-ans[[9]] 

# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(2,2))
############################# tau(q) plot ############################################
plot(qseq,tau, main="Rényi exponent",cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "q")

############################# F(s) plot ##############################################
plot(log10(rs[1:nr]),log10(dmse[1:nr]), main="Fluctuation",
	cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "log(s)", ylab = "logF(s)",
	ylim=c(log10(min(dmse[dmse>0])),log10(max(dmse))))
for (i in 2:nq){ 
  points(log10(rs[1:nr]), log10(dmse[(i*MAX_BOX+1):(i*MAX_BOX+nr)]), col=i,pch=1) 
}
########################## calculate theoretical f(alpha) curve ######################
tot<-100
f_t <- numeric(tot)
alpha_t <- numeric(tot)
for (i in 1:tot){
q<-qmin+i*(qmax-qmin)/tot
f_t[i]<- 0.2000000000e-9 * (0.2075187497e10 * q * exp(-0.2876820725e0 * q) + 
	0.1000000000e11 * q * exp(-0.1386294361e1 * q) + 
	0.7213475205e10 * log(exp(-0.2876820725e0 * q) + 
	exp(-0.1386294361e1 * q)) * exp(-0.2876820725e0 * q) + 
	0.7213475205e10 * log(exp(-0.2876820725e0 * q) + 
	exp(-0.1386294361e1 * q)) * exp(-0.1386294361e1 * q)) / (exp(-0.2876820725e0 * q) + 
	exp(-0.1386294361e1 * q))

alpha_t[i]<- 0.1e1 / q - log((0.75e0^q) + (0.25e0^q)) / q / log(0.2e1) + 
	q * (-(q^-0.2e1) - (-0.2876820725e0 * (0.75e0^q) - 
	0.1386294361e1 * (0.25e0^q)) / ((0.75e0^q) + 
	(0.25e0^q)) / q / log(0.2e1) + log((0.75e0^q) + 
	(0.25e0^q)) * (q^-0.2e1) / log(0.2e1))
}

############################# f(alpha) plot ############################################
plot(alpha,f, main="Multifractal spectrum",cex=1, cex.lab=1.6, cex.axis=1.6, col=4)
######### theoretcal curve in red
lines(alpha_t,f_t,type="l", col=2, lwd=2)

############################# F(q) plot ############################################
plot(qseq,H, main="Generalized Hurst exponent",cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "q")


