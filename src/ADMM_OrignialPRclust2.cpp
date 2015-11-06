#include <Rcpp.h>
using namespace Rcpp;

// calculate distance
// [[Rcpp::export]]
double distance_2(NumericMatrix data, int ndim, int j)
{
    double out  = 0;
    for (int i = 0 ; i<ndim; i++) {
        out = out + pow(data(i,j),2);
    }
    out = sqrt(out);
    return out;
}


// calculate the distance for ADMM
// [[Rcpp::export]]
double distance_umu(NumericMatrix u, NumericMatrix data, int ndim, int i,int j,int uj)
{
    double out  =0;
    for (int ii = 0 ; ii<ndim; ii++) {
        out = out + pow((data(ii,i)-data(ii,j))-u(ii,uj),2);
    }
    out = sqrt(out);
    return out;
}


// calculate distance
// [[Rcpp::export]]
double residual_mu(NumericMatrix mu1, NumericMatrix mu, int ndim,  int numbers)
{
    double numerator  =0;
    double denominator  =0;
    double out  =0;
    for (int i = 0 ; i<ndim; i++) {
        for (int j = 0; j < numbers; j++) {
            numerator = numerator + pow((mu(i,j)-mu1(i,j)),2);
            denominator = denominator + pow(mu1(i,j),2);
        }
    }
    numerator = sqrt(numerator);
    denominator = sqrt(denominator);
    
    out = numerator/denominator;
    return out;
}


// Jude if the theta[i][j] equals 0
// out = -1 theta[i][j] != 0
// out = 1 theta[i][j] = 0
// [[Rcpp::export]]
int is_zero_theta(NumericMatrix theta, int j, int ndim)
{
    int out = -1;
    int ncount = 0;
    for (int i = 0; i<ndim; i++) {
        if (theta(i,j)==0)
        {
            ncount++;
        }
    }
    if (ncount == ndim) {
        out = 1;
    }
    return out;
}


//stopping_criteria(mu,mu1, ndim, numbers, count)
// stopping criteria
// [[Rcpp::export]]
int stopping_criteria(NumericMatrix mu, NumericMatrix mu1, int ndim, int numbers, int count, double epsilon)
{
    int out = 0;
    if (count<2) {
        out = 0;
    } else {
        double Residual;
        double tol = epsilon;
        Residual = residual_mu(mu,mu1,ndim,numbers);
        if (Residual < tol) {
            out = -1;
        }
    }
    
    if (count > 3000) {
        out = -1;
    }
    
    return out;
}
// method =0 PRcluster
// method = 1 Lasso
// rho : the step size
// lambda2: the tuning parameter for theta
// [[Rcpp::export]]
List PRclustADMM(NumericMatrix data, double rho, double lambda2, double tau, int mumethod = 0,  int methods = 0, double epsilon = 0.001)
{
    int count = 0;
    int ndim = data.nrow();
    int numbers = data.ncol();
    int theta_num = numbers*numbers;
    
    
    NumericMatrix mu(ndim,numbers);
    NumericMatrix theta(ndim,theta_num);
    NumericMatrix u(ndim,theta_num);
    NumericMatrix mu1(ndim,numbers);
    NumericMatrix theta1(ndim,theta_num);
    // set the initial value for mu and theta
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<numbers; j++) {
            mu(i,j) = data(i,j);
        }
    }
    
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            theta(i,j) = mu(i,row)-mu(i,col);
        }
    }
    
    
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            u(i,j) = 0;
        }
    }
    double lambda1 = rho;
  
    NumericVector group(numbers);
    while (stopping_criteria(mu,mu1, ndim, numbers, count,epsilon)==0) {
        
        // update mu
        if (mumethod == 0) {
            for (int j=0; j<numbers; j++) {
                for (int i = 0; i<ndim; i++) {
                    double temp_sum1=0;
                    double temp_sum2=0;
                    for (int k = (j+1); k<numbers;k++)
                    {
                        
                        int temp_theta_num = 0;
                        temp_theta_num = j*numbers + k;
                        temp_sum1 = temp_sum1 + mu(i,k) + theta(i,temp_theta_num) + u(i,temp_theta_num);
                    }
                    for (int k = 0; k<j;k++)
                    {
                        
                        int temp_theta_num = 0;
                        temp_theta_num = k*numbers + j;
                        temp_sum2 = temp_sum2 + mu1(i,k) - theta(i,temp_theta_num) - u(i,temp_theta_num);
                    }
                    
                    mu1(i,j) = mu(i,j);
                    
                    double temp_mu =  (data(i,j)+lambda1*temp_sum1+lambda1*temp_sum2) / (1 + lambda1*(numbers-1));
                    mu(i,j) = temp_mu;
                }
            }
        }
        // l_1
        if (mumethod == 1) {
            for (int j=0; j<numbers; j++) {
                for (int i = 0; i<ndim; i++) {
                    double temp_sum1=0;
                    double temp_sum2=0;
                    for (int k = (j+1); k<numbers;k++)
                    {
                        int temp_theta_num = 0;
                        temp_theta_num = j*numbers + k;
                        temp_sum1 = temp_sum1 + mu(i,k) + theta(i,temp_theta_num) + u(i,temp_theta_num) - data(i,j);
                    }
                    for (int k = 0; k<j;k++)
                    {
                        int temp_theta_num = 0;
                        temp_theta_num = k*numbers + j;
                        temp_sum2 = temp_sum2 + mu1(i,k) - theta(i,temp_theta_num) - u(i,temp_theta_num) - data(i,j);
                    }
                    
                    mu1(i,j) = mu(i,j);
                    
                    double temp_lambda;
                    temp_lambda = 1 / (2 * lambda1 * (numbers-1));
                    
                    double temp_mu =  (temp_sum1 + temp_sum2) / (numbers-1);
                    double temp_2 ;
                    int mu_sign = 1;
                    if (temp_mu < 0) {
                        mu_sign = -1;
                    }
                    if (temp_mu == 0){
                        mu_sign = 0;
                    }
                    
                    temp_2 = mu_sign * temp_mu - temp_lambda;
                    if (temp_2 < 0) {
                        temp_2 = 0;
                    }

                    mu(i,j)= mu_sign * temp_2 + data(i,j);
                }
            }
        }
        
        // update theta and  u
        // update theta
        // method == 0, PRcluster
        int theta_num = numbers*numbers;
        if (methods == 0) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    if (distance_2(theta, ndim, j)>= tau) {
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mu(i,row)-mu(i,col)-u(i,j);
                        }
                    } else {
                        double temp,temp1;
                        temp = distance_umu(u,mu, ndim, row, col,j)-lambda2/lambda1;
                        temp1 = distance_umu(u,mu, ndim, row, col,j);
                        if (temp>=0) {
                        } else {
                            temp = 0;
                        }
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= (mu(i,row)-mu(i,col)-u(i,j)) * temp / temp1;
                        }
                    }
                }
            }
        }
        
        // methods== 1 Lasso
        if (methods == 1) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    for (int i = 0; i<ndim; i++) {
                        double abstheta = theta(i,j);
                        if (theta(i,j) <0) {
                            abstheta = -theta(i,j);
                        }
                        if (abstheta >= tau) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mu(i,row)-mu(i,col)-u(i,j);
                        } else {
                            double temp_lambda;
                            temp_lambda = lambda2/lambda1;
                            double lasso_temp ;
                            lasso_temp = mu(i,row)-mu(i,col)-u(i,j);
                            double temp_2 ;
                            int theta_sign = 1;
                            if (lasso_temp < 0) {
                                theta_sign = -1;
                            }
                        
                            if (lasso_temp == 0){
                                theta_sign = 0;
                            }
                        
                            temp_2 = theta_sign*lasso_temp - temp_lambda;
                            if (temp_2 < 0) {
                                temp_2 = 0;
                            }
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= theta_sign * temp_2;
                        }
                    }
                }
            }
        }
        
        // mcp
        if (methods == 2) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                // this may be wrong, I need double check the formula
                if (row <col) {
                    if (distance_2(theta, ndim, j)> (tau * lambda2/lambda1)) {
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mu(i,row)-mu(i,col)-u(i,j);
                        }
                    } else {
                        double temp,temp1;
                        temp = distance_umu(u,mu, ndim, row, col,j)-lambda2/lambda1;
                        temp1 = distance_umu(u,mu, ndim, row, col,j);
                        if (temp>=0) {
                        } else {
                            temp = 0;
                        }
                        double mcp = tau /(tau-1);
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mcp * (mu(i,row)-mu(i,col)-u(i,j)) * (temp / temp1);
                        }
                    }
                }
            }
        }
        
        // grouped SCAD
        if (methods == 3) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    if (distance_2(theta, ndim, j)> (lambda2 * tau/lambda1)) {
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mu(i,row)-mu(i,col)-u(i,j);
                        }
                    } else if (distance_2(theta, ndim, j) <= (lambda2 * tau/lambda1) && distance_2(theta, ndim, j) > (2*lambda2/lambda1) ) {
                        
                        double temp,temp1;
                        temp = distance_umu(u,mu, ndim, row, col,j)- (tau/(tau-1)) * lambda2/lambda1;
                        temp1 = distance_umu(u,mu, ndim, row, col,j);
                        if (temp <0) {
                            temp = 0;
                        }
                        double scad = (tau -1)/(tau-2);
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= scad* (mu(i,row)-mu(i,col)-u(i,j)) * (temp / temp1);
                        }
                        
                    } else {
                        double temp,temp1;
                        temp = distance_umu(u,mu, ndim, row, col,j)-lambda2/lambda1;
                        temp1 = distance_umu(u,mu, ndim, row, col,j);
                        if (temp>=0) {
                        } else {
                            temp = 0;
                        }
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= (mu(i,row)-mu(i,col)-u(i,j)) * temp / temp1;
                        }
                    }
                }
            }
        }
        
        // update u
        for (int j=0; j<theta_num; j++) {
            int row = j / numbers;
            int col = j % numbers;
            // only update the theta we need.
            if (row <col) {
                for(int i = 0;i <ndim;i++)
                {
                    u(i,j) =u(i,j) + theta(i,j)-(mu(i,row)-mu(i,col));
                }
            }
        }

        count++;     
    }
    
    // get the final group
    int group_num = 1;
    for (int i = 0; i < numbers; i++) {
        group[i] = -1;
    }
    group[0] =group_num;
    for (int i = 0; i<numbers; i++) {
        int temp_start = i+1;
        for (int j = temp_start; j<numbers; j++) {
            int temp = i*numbers + j;
            // use the following function to determine whether theta[i][j] equals zero or not;
            if (is_zero_theta(theta, temp, ndim) == 1) {
                if (group[j] == -1) {
                    if (group[i] != -1) {
                        group[j] = group[i];
                    } else {
                        group_num++;
                        group[i] = group_num;
                        group[j] = group[i];
                    }
                } else {
                    if (group[i] == -1) {
                        group[i] = group[j];
                    } else {
                        if (group[i] < group[j]) {
                            for (int tem_i = 0; tem_i < numbers; tem_i ++) {
                                if (group[tem_i] == group[j]) {
                                    group[tem_i] = group[i];
                                }
                            }
                        } else if (group[i] > group[j]){
                            for (int tem_i = 0; tem_i < numbers; tem_i ++) {
                                if (group[tem_i] == group[i]) {
                                    group[tem_i] = group[j];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < numbers; i++) {
        if (group[i] == -1) {
            group_num++;
            group[i] = group_num;
        }
    }
    List out;
    out["mu"] =mu;
    out["theta"] = theta;
    out["data"] = data;
    out["count"] = count;
    out["group"] = group;
    return out;
}


// [[Rcpp::export]]
List clusterStat(NumericVector trueGroup,NumericVector group)
{
    // return rand, adjusted rand and Jaccard
    int numbers = trueGroup.size();
    
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    
    for (int i = 0; i < (numbers- 1); i ++) {
        for (int j = (i + 1); j < numbers; j ++) {
            if (trueGroup[i] == trueGroup[j] &&  group[i] == group[j]) {
                a = a + 1;
            }
            if (trueGroup[i] != trueGroup[j] &&  group[i] != group[j]) {
                b = b + 1;
            }
            if (trueGroup[i] == trueGroup[j] &&  group[i] != group[j]) {
                c = c + 1;
            }
            if (trueGroup[i] != trueGroup[j] &&  group[i] == group[j]) {
                d = d + 1;
            }
            
        }
    }
    
    //Jaccard = a/( a+c+d);
    List out;
    //out["Rand"]= rand;
    //out["Jaccard"] =Jaccard;
    //out["number"] = numbers;
    out["a"] =a;
    out["b"] =b;
    out["c"] =c;
    out["d"] =d;
    return out;
}



// Original Method

// calculate distance
// [[Rcpp::export]]
double distance_mu(NumericMatrix data, int ndim, int i,int j)
{
    double out  =0;
    for (int ii = 0 ; ii<ndim; ii++) {
        out = out + pow((data(ii,i)-data(ii,j)),2);
    }
    out = sqrt(out);
    return out;
}


// calculate the S
// [[Rcpp::export]]
double cal_S(NumericMatrix data,NumericMatrix mu, NumericMatrix theta, double lambda1, double lambda2, double tau,int ndim, int numbers, int methods)
{
    double out;
    if (methods == 0) {
        double S1;
        double S2;
        double temp1=0, temp2=0, temp3=0, temp4=0;
        for (int j = 0; j<numbers; j++) {
            double dist=0;
            for (int i = 0; i<ndim; i++) {
                dist = dist + pow((data(i,j)-mu(i,j)),2);
            }
            temp1 = temp1 + dist;
        }
        int theta_num = numbers*numbers;
        
        for (int j = 0 ; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            // only sum the theta we need.
            
            if (row <col) {
                double dist1 = 0;
                double dist2 = 0;
                double dist3 = 0;
                for (int i =0; i<ndim; i++) {
                    dist1 = dist1 + pow((mu(i,row)-mu(i,col)-theta(i,j)), 2);
                    dist2 = dist2 + pow(theta(i,j), 2);
                }
                temp2 = temp2 + dist1;
                temp3 = temp3 + sqrt(dist2);
                dist3 = sqrt(dist2)-tau;
                if (dist3 >= 0) {
                    temp4 = temp4 + dist3;
                }
            }
        }
        S1 = 0.5 * temp1 + lambda1* 0.5 * temp2 + lambda2 * temp3;
        S2 =  lambda2 * temp4;
        out = S1 - S2;
    }
    if (methods == 1) {
        double S1;
        double temp1=0, temp2=0, temp3=0;
        for (int j = 0; j<numbers; j++) {
            double dist=0;
            for (int i = 0; i<ndim; i++) {
                dist = dist + pow((data(i,j)-mu(i,j)),2);
            }
            temp1 = temp1 + dist;
        }
        int theta_num = numbers*numbers;
        
        for (int j = 0 ; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            // only sum the theta we need.
            
            if (row < col) {
                double dist1 = 0;
                double dist2 = 0;
                for (int i = 0; i<ndim; i++) {
                    dist1 = dist1 + pow((mu(i,row)-mu(i,col)-theta(i,j)), 2);
                    if (theta(i,j)>0) {
                        dist2 = dist2 + theta(i,j);
                    }
                    else
                    {
                        dist2 = dist2 - theta(i,j);
                    }
                    
                }
                temp2 = temp2 + dist1;
                temp3 = temp3 + dist2;
            }
        }
        
        S1 = 0.5 * temp1 + lambda1* 0.5 * temp2 + lambda2 * temp3;
        out = S1;
    }
    return out;
}



// determine whether terminate iteration or not
// out = 0 : continue;
// out = -1 : terminate.
// also need change the judge function since different methods use different S;
// method =0 PRcluster
// method = 1 Lasso
// [[Rcpp::export]]
int judge_iteration(NumericMatrix data,NumericMatrix mu, NumericMatrix theta,NumericMatrix mu1, NumericMatrix theta1, double lambda1, double lambda2, double tau,int ndim, int numbers, int count, int methods)
{
    int out = 0;
    double S1;
    double S2;
    if (count<1) {
        out = 0;
    }
    else{
        S1 = cal_S(data,mu,theta,lambda1,lambda2,tau,ndim,numbers,methods);
        S2 = cal_S(data,mu1,theta1,lambda1,lambda2,tau,ndim,numbers,methods);
        if ((S1 - S2)>0 ) {
            out = -1;
        }
    }
    return out;
}


// method =0 PRcluster
// method = 1 Lasso
// [[Rcpp::export]]
List PRclustOriginal(NumericMatrix data, double lambda1, double lambda2, double tau, int mumethod = 0, int methods = 0)
{
    int count = 0;
    int ndim = data.nrow();
    int numbers = data.ncol();
    int theta_num = numbers*numbers;
    
    
    NumericMatrix mu(ndim,numbers);
    NumericMatrix theta(ndim,theta_num);
    
    // set the initial value for mu and theta
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<numbers; j++) {
            mu(i,j) = data(i,j);
        }
    }
    
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            theta(i,j) = mu(i,row)-mu(i,col);
        }
    }
    
    
    // define mu1 theta1 and group
    NumericMatrix mu1(ndim,numbers);
    NumericMatrix theta1(ndim,theta_num);
    NumericVector group(numbers);
    //judge_iteration(data,mu,theta,mu1,theta1,lambda1,lambda2,tau,ndim,numbers,count,methods) ==0
    
    while (judge_iteration(data,mu,theta,mu1,theta1,lambda1,lambda2,tau,ndim,numbers,count,methods) ==0) {
        // update mu
        for (int j=0; j<numbers; j++) {
            for (int i = 0; i<ndim; i++) {
                double temp_sum1=0;
                double temp_sum2=0;
                for (int k = (j+1); k<numbers;k++)
                {
                    
                    int temp_theta_num = 0;
                    temp_theta_num = j*numbers + k;
                    temp_sum1 = temp_sum1 + mu(i,k) + theta(i,temp_theta_num);
                }
                for (int k = 0; k<j;k++)
                {
                    
                    int temp_theta_num = 0;
                    temp_theta_num = k*numbers + j;
                    temp_sum2 = temp_sum2 + mu(i,k) - theta(i,temp_theta_num);
                }
                mu1(i,j) = mu(i,j);
                mu(i,j) = (data(i,j)+lambda1*temp_sum1+lambda1*temp_sum2)/(1+lambda1*(numbers-1));
            }
        }
        
        // update theta
        // method == 0, PRcluster
        if (methods == 0) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    if (distance_2(theta, ndim, j)>= tau) {
                        for (int i = 0; i<ndim; i++) {
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= mu(i,row)-mu(i,col);
                        }
                    } else {
                        for (int i = 0; i<ndim; i++) {
                            double temp,temp1;
                            temp = distance_mu(mu, ndim, row, col)-lambda2/lambda1;
                            temp1 = distance_mu(mu, ndim, row, col);
                            if (temp>=0) {
                            } else {
                                temp = 0;
                            }
                            theta1(i,j) = theta(i,j);
                            theta(i,j)= (mu(i,row)-mu(i,col))*temp/temp1;
                        }
                    }
                }
            }
        }
        
        // methods== 1 Lasso
        if (methods == 1) {
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    for (int i = 0; i<ndim; i++) {
                        double temp_lambda;
                        temp_lambda = lambda2/lambda1;
                        double lasso_temp ;
                        lasso_temp = mu(i,row)-mu(i,col);
                        double temp_2 ;
                        int theta_sign = 1;
                        if (lasso_temp < 0) {
                            theta_sign = -1;
                        }
                        
                        if (lasso_temp == 0){
                            theta_sign = 0;
                        }
                        
                        temp_2 = theta_sign*lasso_temp - temp_lambda;
                        if (temp_2 < 0) {
                            temp_2 = 0;
                        }
                        theta1(i,j) = theta(i,j);
                        theta(i,j)= theta_sign * temp_2;
                    }
                }
            }
        }
        count++;
    }
    
    // get the final group
    int group_num = 1;
    for (int i = 0; i < numbers; i++) {
        group[i] = -1;
    }
    group[0] =group_num;
    for (int i = 0; i<numbers; i++) {
        int temp_start = i+1;
        for (int j = temp_start; j<numbers; j++) {
            int temp = i*numbers + j;
            // use the following function to determine whether theta[i][j] equals zero or not;
            if (is_zero_theta(theta, temp, ndim) == 1) {
                if (group[j] == -1) {
                    if (group[i] != -1) {
                        group[j] = group[i];
                    } else {
                        group_num++;
                        group[i] = group_num;
                        group[j] = group[i];
                    }
                } else {
                    if (group[i] == -1) {
                        group[i] = group[j];
                    } else {
                        if (group[i] < group[j]) {
                            for (int tem_i = 0; tem_i < numbers; tem_i ++) {
                                if (group[tem_i] == group[j]) {
                                    group[tem_i] = group[i];
                                }
                            }
                        } else if (group[i] > group[j]){
                            for (int tem_i = 0; tem_i < numbers; tem_i ++) {
                                if (group[tem_i] == group[i]) {
                                    group[tem_i] = group[j];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < numbers; i++) {
        if (group[i] == -1) {
            group_num++;
            group[i] = group_num;
        }
    }
    List out;
    out["mu"] =mu;
    out["theta"] = theta;
    out["data"] = data;
    out["count"] = count;
    out["group"] = group;
    return out;
}


