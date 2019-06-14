#data step
data_NPI<-read.xlsx("NPI.xlsx",sheetName='Sheet1', colIndex=1)
data_FFrmrf<-read.table('FF_rmrf.txt')
data_RF<-read.table('FF_rf.txt')

#-----------question 3.1-------------
NPI_ts <- ts(data_NPI, start=c(1978, 1), end=c(2015, 1), frequency=4)

plot(NPI_ts, main='Plot for NPI data', xlab='time', ylab='NPI')
#classic decomposition of the NPI
NPI_ts_decompose<- decompose(NPI_ts)
plot(NPI_ts_decompose)
mean(NPI_ts)
var(NPI_ts)
par(mfrow=c(2,1))
acf(NPI_ts, main='ACF plot for NPI data', lag.max=149)
pacf(NPI_ts, main='PACF plot for NPI data', lag.max=149)
sharp_ratio_smooth<-c()
for(i in 21:149){
  sharp_ratio_i<-mean(data_NPI[i,]-data_RF[i,])/sd(data_NPI[1:i,])
  sharp_ratio_smooth[i-20]<-sharp_ratio_i
}
par(mfrow=c(1,1))
plot(sharp_ratio_smooth, main='sharp ratio for smooth NPI', ylab='sharp ratio')
summary(sharp_ratio_smooth)
FFrmrf_ts <- ts(data_FFrmrf, start=c(1978, 1), end=c(2015, 1), frequency=4)
plot(FFrmrf_ts, main='Plot for FFrm-rf data', xlab='time', ylab='FFrm-rf')
FFrmrf_ts_decompose<- decompose(FFrmrf_ts)
plot(FFrmrf_ts_decompose)
par(mfrow=c(2,1))
acf(FFrmrf_ts, main='ACF plot for FFrm-rf data', lag.max=50)
pacf(FFrmrf_ts, main='PACF plot for FFrm-rf data', lag.max=50)
mean(FFrmrf_ts)
var(FFrmrf_ts)
#--------------question3.2--------------
#MA(q) process
AIC_MA<-c()
for(i in 1:15){
  AIC_MA[i]<-arima(NPI_ts, c(0,0,i))$aic}
AIC_MA_min<-min(AIC_MA)
MA_order<-which.min(AIC_MA)
MA_model<-arima(NPI_ts-mean(NPI_ts), c(0,0,MA_order))

#conduct LB and BP test for MA model
BoxPierce(MA_model, lags=c(10, 15, 20))
LjungBox(MA_model, lags=c(10, 15, 20))

#retrieve r_t
MAmodel_coef<-MA_model$coef[1:8]
c<-1+sum(MAmodel_coef)
MAr_t<- c*MA_model$residuals + mean(NPI_ts)
par(mfrow=c(1,1))
plot(MAr_t, main='Plot for unsmooth NPI (using MA(8) process)', ylab='r_t')
summary(MAr_t)
par(mfrow=c(2,1))
acf(MAr_t, main='ACF plot for unsmoothed NPI data using MA(8) process', lag.max=50)
pacf(MAr_t, main='PACF plot for unsmoothed NPI data using MA(8) process', lag.max=50)
#volatility
sd(MAr_t)

#AR(1) process
AR_1_model<-arima(NPI_ts-mean(NPI_ts), c(1,0,0))
w0_hat_1<-1-sum(AR_1_model$coef[1])
ARr_t_1<-AR_1_model$residuals/w0_hat_1 + mean(NPI_ts)
plot(ARr_t_1, main='Plot for unsmoothed NPI (using AR(1) process)', ylab='r_t')
summary(ARr_t_1)
par(mfrow=c(2,1))
acf(ARr_t_1, main='ACF plot for unsmoothed NPI data using AR(1) process', lag.max=50)
pacf(ARr_t_1, main='PACF plot for unsmoothed NPI data using AR(1) process', lag.max=50)
#conduct LB and BP test for AR1 model
BoxPierce(AR_1_model, lags=c(10, 15, 20))
LjungBox(AR_1_model, lags=c(10, 15, 20))

#AR(1,4) process
fit<-lm(NPI_ts-mean(NPI_ts) ~ lag(NPI_ts-mean(NPI_ts),-1)+lag(NPI_ts-mean(NPI_ts),-4))
AR14_model<-dyn$lm(NPI_ts-mean(NPI_ts) ~ lag(NPI_ts-mean(NPI_ts),-1)+lag(NPI_ts-mean(NPI_ts),-4))
w0_hat_14<-1-sum(AR14_model$coef[2:3])
ARr_t_14<-AR14_model$residuals/w0_hat_14 + mean(NPI_ts)
plot(ARr_t_14, main='Plot for unsmoothed NPI (using AR(1, 4) process)', ylab='r_t', type='l')
summary(ARr_t_14)
par(mfrow=c(2,1))
acf(ARr_t_14, main='ACF plot for unsmoothed NPI data using AR(1,4) process', lag.max=50)
pacf(ARr_t_14, main='PACF plot for unsmoothed NPI data using AR(1,4) process', lag.max=50)
#conduct LB and BP test for AR full model
Box.test(AR14_model$residuals, lag=10, type = "Ljung-Box", fitdf = 0)
Box.test(AR14_model$residuals, lag=15, type = "Ljung-Box", fitdf = 0)
Box.test(AR14_model$residuals, lag=20, type = "Ljung-Box", fitdf = 0)
Box.test(AR14_model$residuals, lag=10, type = "Box-Pierce", fitdf = 0)
Box.test(AR14_model$residuals, lag=15, type = "Box-Pierce", fitdf = 0)
Box.test(AR14_model$residuals, lag=20, type = "Box-Pierce", fitdf = 0)

#AR(p) process
AIC_AR<-c()
for(j in 1:15){
  AIC_AR[j]<-arima(NPI_ts, c(j,0,0))$aic
}
AIC_AR_min<-min(AIC_AR)
AR_order<-which.min(AIC_AR)
AR_full_model<-arima(NPI_ts-mean(NPI_ts), c(AR_order,0,0))

#conduct LB and BP test for AR full model
BoxPierce(AR_full_model, lags=c(10, 15, 20))
LjungBox(AR_full_model, lags=c(10, 15, 20))

#retrieve r_t using AR full
w0_hat_p<-1-sum(AR_full_model$coef[1:7])
ARr_t_full<-AR_full_model$residuals/w0_hat_p +mean(NPI_ts)
plot(ARr_t_full, main='Plot for unsmooth NPI (using AR(7) process)', ylab='r_t')
summary(ARr_t_full)
par(mfrow=c(2,1))
acf(ARr_t_full, main='ACF plot for unsmoothed NPI data using AR(7) process', lag.max=50)
pacf(ARr_t_full, main='PACF plot for unsmoothed NPI data using AR(7) process', lag.max=50)
par(mfrow=c(1,1))
#sharp ratio for MA, unsmooth
sharp_ratio_unsmooth_ma<-c()
for(i in 21:149){
  sharp_ratio_ma_i<-mean(MAr_t[i]-data_RF[i,])/sd(MAr_t[1:i])
  sharp_ratio_unsmooth_ma[i-20]<-sharp_ratio_ma_i
}

plot(sharp_ratio_unsmooth_ma, main='sharpe ratio for unsmooth NPI using MA process', ylab='sharp ratio')

#sharp ratio for AR full, unsmooth
sharp_ratio_unsmooth_ar<-c()
for(i in 21:149){
  sharp_ratio_ar_i<-mean(ARr_t_full[i]-data_RF[i,])/sd(ARr_t_full[1:i])
  sharp_ratio_unsmooth_ar[i-20]<-sharp_ratio_ar_i
}
plot(sharp_ratio_unsmooth_ar, main='sharpe ratio for unsmooth NPI using AR(7) process', ylab='sharp ratio')

#sharp ratio for AR14, unsmooth
sharp_ratio_unsmooth_ar14<-c()
for(i in 21:145){
  sharp_ratio_ar14_i<-mean(ARr_t_14[i]-data_RF[i,])/sd(ARr_t_14[1:i])
  sharp_ratio_unsmooth_ar14[i-20]<-sharp_ratio_ar14_i
}
plot(sharp_ratio_unsmooth_ar14, main='sharpe ratio for unsmooth NPI using AR(1,4) process', ylab='sharp ratio')
#sharp ratio for AR14, unsmooth
sharp_ratio_unsmooth_ar1<-c()
for(i in 21:149){
  sharp_ratio_ar1_i<-mean(ARr_t_1[i]-data_RF[i,])/sd(ARr_t_1[1:i])
  sharp_ratio_unsmooth_ar1[i-20]<-sharp_ratio_ar1_i
}
plot(sharp_ratio_unsmooth_ar1, main='sharpe ratio for unsmooth NPI using AR(1) process', ylab='sharp ratio')

#-----------------(= 3 =)------------------------
data_SMB<-read.table('FF_SMB.txt')
data_HML<-read.table('FF_HML.txt')
RF_ts <- ts(data_RF, start=c(1978, 1), end=c(2015, 1), frequency=4)
SMB_ts <- ts(data_SMB, start=c(1978, 1), end=c(2015, 1), frequency=4)
HML_ts <- ts(data_HML, start=c(1978, 1), end=c(2015, 1), frequency=4)

#using NPI data
FLM_NPI<-lm(NPI_ts-RF_ts ~ SMB_ts+HML_ts+FFrmrf_ts) #FLM short for factor loading model
Box.test(FLM_NPI$residuals, lag=10, type = "Ljung-Box", fitdf = 0)
#using MA(q) 
FLM_MA<-lm(MAr_t-RF_ts ~ SMB_ts+HML_ts+FFrmrf_ts)
plot(FLM_MA)
Box.test(FLM_MA$residuals, lag=10, type = "Ljung-Box", fitdf = 0)

#using AR(1,4)
FLM_ARr14<-lm(ARr_t_14-RF_ts[1:145] ~ SMB_ts[1:145]+HML_ts[1:145]+FFrmrf_ts[1:145])
plot(FLM_ARr14)
Box.test(FLM_ARr14$residuals, lag=10, type = "Ljung-Box", fitdf = 0)
#using AR(P)
FLM_ARrp<-lm(ARr_t_full-RF_ts ~ SMB_ts+HML_ts+FFrmrf_ts)
plot(FLM_ARrp)
Box.test(FLM_ARrp$residuals, lag=10, type = "Ljung-Box", fitdf = 0)
#using AR(1)
FLM_ARr1<-lm(ARr_t_1-RF_ts ~ SMB_ts+HML_ts+FFrmrf_ts)
plot(FLM_ARr1)
#PICMO method
omega<-MAmodel_coef/c

#SMB transformed factors
smb_fct<-c()
for(i in 9:149){
  product<-c()
  for(j in 1:9){if (j == 1) product[j]<-(1/c)*SMB_ts[i-j+1]
                else product[j]<-omega[j-1]*SMB_ts[i-j+1]}
  smb_fct[i-8]<-sum(product)}

#HMLtransformed factors
hml_fct<-c()
for(i in 9:149){
  product<-c()
  for(j in 1:9){
    if (j == 1) product[j]<-(1/c)*HML_ts[i-j+1]
    else product[j]<-omega[j-1]*HML_ts[i-j+1]}
  hml_fct[i-8]<-sum(product)}

#Mkt transformed factors
mkt_fct<-c()
for(i in 9:149){
  product<-c()
  for(j in 1:9){if (j == 1) product[j]<-(1/c)*FFrmrf_ts[i-j+1]
                else product[j]<-omega[j-1]*FFrmrf_ts[i-j+1]}
                mkt_fct[i-8]<-sum(product)}

PICMO<-lm(NPI_ts[9:149]~hml_fct+mkt_fct+smb_fct)
PICMO
Box.test(PICMO$residuals, lag=10, type = "Ljung-Box", fitdf = 0)
