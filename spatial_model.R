lamda = 2;
c = 1;
a = 0.068;
dh = 0.55;
dp = 0.25;

T = 120 #how many itarations
n = 20 #size of grid
vectorH = rep(0,T);
vectorP = rep(0,T);
vector1 = rep(0,n);
vector2 = rep(0,n);
H = array(c(vector1, vector2), dim = c(n, n, 1)); 
P = array(c(vector1, vector2), dim = c(n, n, 1));
Hpast = array(c(vector1, vector2), dim = c(n, n, 1)); 
Ppast = array(c(vector1, vector2), dim = c(n, n, 1)); 
Haux = array(c(vector1, vector2), dim = c(n, n, 1)); 
Paux = array(c(vector1, vector2), dim = c(n, n, 1)); 

mysum <- function(size,i,j,arr){
  result=0;
  print(j);
  if(i == 1 && j == 1){ #top left corner
    result = arr[i,j+1,1] + arr[i+1,j+1,1] + arr[i+1,j,1] + arr[i+1,size,1] + arr[i,size,1] + arr[size,size,1] + arr[size,j,1] + arr[size,j+1,1];
  }
  else if(i == 1 && j ==size){ #top right corner
    result = arr[1,1,1] +   arr[i+1,1,1] +   arr[i+1,j,1] + arr[i+1,j-1,1] + arr[i,j-1,1] +  arr[size,j-1,1] + arr[size,size,1] + arr[size,1,1];
  }
  else if(i==size && j==size){ #bottom right corner
    result = arr[i,1,1] + arr[1,1,1] + arr[1,size,1] + arr[1,size-1,1] + arr[i,j-1,1] + arr[i-1,j-1,1] + arr[i-1,j,1] + arr[i-1,1,1];
  }
  else if(i==size && j == 1){ #bottom left corner
    result = arr[i,j+1,1] + arr[1,j+1,1] + arr[1,1,1] + arr[1,size,1] + arr[size,size,1] + arr[i-1,size,1] + arr[i-1,j,1] + arr[i-1,j+1,1];
  }
  else if( i == 1){ # upper bound
    result = arr[i,j+1,1] + arr[i+1,j+1,1] + arr[i+1,j,1] + arr[i+1,j-1,1] + arr[i,j-1,1] + arr[size,j-1,1] +arr[size,j,1] + arr[size,j+1,1];
  }
  else if(i == size){ # lower bound
    result = arr[i,j+1,1] + arr[1,j+1,1] + arr[1,j,1] + arr[1,j-1,1] +arr[i,j-1,1] + arr[i-1,j-1,1] + arr[i-1,j,1] +arr[i-1,j+1,1];
  }
  else if(j == 1){ # left bound
    result = arr[i,j+1,1] + arr[i+1,j+1,1] +arr[i+1,j,1] + arr[i+1,size,1] + arr[i,size,1] + arr[i-1,size,1] + arr[i-1,j,1] + arr[i-1,j+1,1];
  }
  else if(j == size){ # right bound
    result = arr[i,1,1] + arr[i+1,1,1] + arr[i+1,j,1] + arr[i+1,j-1,1] + arr[i,j-1,1] + arr[i-1,j-1,1] + arr[i-1,j,1] + arr[i-1,1,1];
  }
  else{ #calc inner cell
    result = arr[i,j+1,1] + arr[i+1,j+1,1] + arr[i+1,j,1] + arr[i+1,j-1,1] + arr[i,j-1,1] + arr[i+1,j-1,1] + arr[i-1,j,1] + arr[i-1,j+1,1];
  }
  return(result);
}
#################### Initilize POPULATION
  for (row in 1:n){
    for (col in 1:n){
      H[row,col,1] = 0;
      Hpast[row,col,1] = 0;
    }
  }
H[n/2,n/2,1] = lamda*log(lamda)/((lamda-1)*a*c);
Hpast[n/2,n/2,1] = lamda*log(lamda)/((lamda-1)*a*c);

P[n/2,n/2,1] = log(lamda)/a;
Ppast[n/2,n/2,1] = log(lamda)/a;
#################### Initilize POPULATION

############## Generate Iterations
for(iteration in 1:T){
  for (row in 1:n){
    for (col in 1:n){
      invisible( capture.output(Haux[row,col,1] <- (1-dh) * Hpast[row,col,1] + dh/8*(mysum(n,row,col,Hpast))));
      invisible(capture.output(Paux[row,col,1] <- (1-dp) * Ppast[row,col,1] + dp/8*(mysum(n,row,col,Ppast))));
      H[row,col,1] = lamda*Haux[row,col,1]*exp(-a*Paux[row,col,1]);
      P[row,col,1] = c*Haux[row,col,1]*(1-exp(-a*Paux[row,col,1]));
      vectorH[iteration] =  vectorH[iteration] + H[row,col,1]; ### Calculate overall population for hosts
      vectorP[iteration] =  vectorP[iteration] + P[row,col,1]; ### Calculate overall population for parasites
    }
  }
  for(row in 1:n){
    for(col in 1:n){
      Hpast[row,col,1] = H[row,col,1];
      Ppast[row,col,1] = P[row,col,1];
    }
  }
}
################ Generate Iterations




mat <- matrix(H, nrow = n, ncol = n)

image(mat,col=grey.colors(n, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL, rev = FALSE))



plot((1:T),vectorP);   #Population for every iteration
plot((1:T),vectorH);
#vectorP
########### Задача 1
    # СКРИПТ    
########### Задача 2
    # ФОРМУЛА ЗА СУМИТЕ ПО ГРАНИЦАТА НА РЕШЕТКАТА
    # 
    #
#dh,dp? за да бъдат популациите > 0
########### Задача 3
    #при dh=0.01 и dp=0.01 се получава спирала

########### Задача 4

########### Задача 5
#  (1)При голямо 'lamda' популацията  от гостоприемници достига по-рано големи популации и се срива след това
#  (2)При голямо 'c' популацията от паразити достига рано пик и след това се срива  
#  (3)При малко 'а' и двете популации достигат пик почти в едни и същи итерации и след това се сриват
#




