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
    
    if (count > 1000) {
        out = -1;
    }
    
    return out;
}


// calculate the intermidate values: ||Ax||_2
// [[Rcpp::export]]
double cal_primalInterValue(NumericMatrix mu, int ndim, int numbers)
{
    double out =0;
    int theta_num = numbers *numbers;
    // the primal resiudal is theta_ij-(mu_i - mu_j), 1<= i<j <N, p= 1,2,...,p
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            if (row <col){
                out = out +pow(mu(i,row)-mu(i,col),2);
                //u[i][j] = 0;
            }
        }
    }
    out = sqrt(out);
    return out;
}

// calculate the intermidate values: ||Ax||_2
// [[Rcpp::export]]
double cal_primalInterValue2(NumericMatrix theta, int ndim, int numbers)
{
    double out =0;
    int theta_num = numbers *numbers;
    // the primal resiudal is theta_ij-(mu_i - mu_j), 1<= i<j <N, p= 1,2,...,p
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            if (row <col){
                out = out +pow(theta(i,j),2);
            }
            //u[i][j] = 0;
        }
    }
    out = sqrt(out);
    return out;
}


// calculate the primal resiudal
// [[Rcpp::export]]
double cal_relResInterValue(NumericMatrix u,double rho, int ndim, int numbers)
{
    double out =0;
    for (int j =0; j<numbers; j++) {
        for (int i = 0; i<ndim; i++) {
            //calculate j and
            //int index;
            double tempout;
            
            if (j == 0) {
                tempout = 0;
                for (int ii = 0; ii < (numbers-1); ii++) {
                    tempout = tempout + u(i,ii+1);
                }
                out = out + pow(tempout, 2);
            }
            else{
                tempout =0;
                for (int ii =0; ii<(numbers-1); ii++) {
                    if (ii < j) {
                        tempout = tempout + u(i,ii*numbers + j);
                    }
                    else{
                        tempout = tempout - u(i,j*numbers + ii +1) ;
                    }
                }
                out = out + pow(tempout,2);
            }
        }
    }
    
    out = rho * sqrt(out);
    return out;
}


// calculate the primal resiudal
// [[Rcpp::export]]
double cal_primalRes(NumericMatrix mu, NumericMatrix theta, int ndim, int numbers)
{
    double out =0;
    int theta_num = numbers *numbers;
    // the primal resiudal is theta_ij-(mu_i - mu_j), 1<= i<j <N, p= 1,2,...,p
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            if (row <col){
                out = out +pow(theta(i,j) - (mu(i,row)-mu(i,col)),2);
            }
            //u[i][j] = 0;
        }
    }
    out = sqrt(out);
    return out;
}

// calcualte the dual resduial
// [[Rcpp::export]]
double cal_dualRes(NumericMatrix theta,NumericMatrix theta1, double rho, int ndim, int numbers)
{
    double out =0;
    for (int j =0; j<numbers; j++) {
        for (int i = 0; i<ndim; i++) {
            //calculate j and
            //int index;
            double tempout;
            
            if (j == 0) {
                tempout = 0;
                for (int ii = 0; ii < (numbers-1); ii++) {
                    tempout = tempout + theta(i,ii+1)- theta1(i,ii+1);
                }
                out = out + pow(tempout, 2);
            }
            else{
                tempout =0;
                for (int ii =0; ii<(numbers-1); ii++) {
                    if (ii < j) {
                        tempout = tempout + theta(i,ii*numbers + j) - theta1(i,ii*numbers + j);
                    }
                    else{
                        tempout = tempout -(theta(i,j*numbers + ii +1)  - theta1(i,j*numbers + ii +1));
                    }
                }
                out = out + pow(tempout,2);
            }
        }
    }
    out =rho * sqrt( out);
    return out;
}



// stopping criteria
// [[Rcpp::export]]
int stopping_criteria2(NumericMatrix mu, NumericMatrix theta, NumericMatrix theta1,NumericMatrix u, double rho, double abs_res, double rel_res, int ndim, int numbers, int count)
{
    int out = 0;
    if (count<1) {
        out = 0;
    }
    else{
        double primalResidual,dualResidual;
        primalResidual = cal_primalRes(mu,theta,ndim,numbers);
        dualResidual = cal_dualRes(theta,theta1,rho,ndim,numbers);
        
        // calculate the Ax;
        double Ax;
        double Bz;
        double isPrimalRes, isDualRes;
        Ax = cal_primalInterValue(mu, ndim, numbers);
        Bz = cal_primalInterValue2(theta, ndim, numbers);
        
        if (Ax > Bz) {
            isPrimalRes = Ax * rel_res + sqrt(ndim*numbers*(numbers-1)/2) * abs_res;
        }
        else{
            isPrimalRes = Bz * rel_res + sqrt(ndim*numbers*(numbers-1)/2) * abs_res;
        }
        
        isDualRes = sqrt(ndim*numbers) * abs_res + rel_res * cal_relResInterValue(u,rho, ndim, numbers);
        // && dualResidual < isDualRes
        if (primalResidual < isPrimalRes && dualResidual < isDualRes ) {
            out = -1;
        }
        
    }
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
    } else {
        S1 = cal_S(data,mu,theta,lambda1,lambda2,tau,ndim,numbers,methods);
        S2 = cal_S(data,mu1,theta1,lambda1,lambda2,tau,ndim,numbers,methods);
        if ((S1 - S2)>0 ) {
            out = -1;
        }
    }
    if (count > 3000) {
        out = -1;
    }
    return out;
}

// For PRclust2
// [[Rcpp::export]]
int judge_iteration2(NumericMatrix PRmu, NumericMatrix mu,double epsilon, int ndim, int numbers,int count)
{
    int out = 0;
    double sum = 0 ;
    double sum2 = 0 ;
    if (count > 1) {
        out = -1;
    }
    else{
        for(int i = 0;i<ndim;i++){
            for (int j =0; j<numbers; j++) {
                sum = sum +  sqrt((PRmu(i,j)-mu(i,j)) * (PRmu(i,j)-mu(i,j)));
                sum2 = sum2 +  sqrt(mu(i,j) * mu(i,j));
            }
        }
        if ((sum/sum2) < epsilon) {
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


// calculate the S
// [[Rcpp::export]]
double cal_S_ADMM(NumericMatrix data, NumericMatrix theta,NumericMatrix theta2, NumericMatrix mu, double lambda2, double tau,int ndim, int numbers)
{
    double S1;
    double temp1=0, temp4=0;
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
            double dist2 = 0;
            double dist3 = 0;
            for (int i =0; i<ndim; i++) {
                dist2 = dist2 + pow(theta2(i,j), 2);
            }
            dist3 = sqrt(dist2)-tau;
            if (dist3 >= 0) {
                temp4 = temp4 + lambda2 * tau;
            } else {
                double dist4 = 0;
                for (int i =0; i<ndim; i++) {
                    dist4 = dist4 + pow(theta(i,j), 2);
                }
                temp4 = temp4 + lambda2 * sqrt(dist4);
            }
        }
    }
    S1 = 0.5 * temp1 + temp4;
    return S1;
}




// rho : the step size
// lambda2: the tuning parameter for theta
// [[Rcpp::export]]
List DCADMM(NumericMatrix data, double rho, double lambda2, double tau ,double abs_res = 0.5, double rel_res = 0.5,int mumethod = 0, int methods =0)
{
    int count = 0;
    int ndim = data.nrow();
    int numbers = data.ncol();
    int theta_num = numbers*numbers;
    
    NumericMatrix mu(ndim,numbers);
    NumericMatrix theta(ndim,theta_num);
    NumericMatrix u(ndim,theta_num);
    NumericMatrix mu1(ndim,numbers);
    NumericMatrix mu2(ndim,numbers);

    NumericMatrix theta1(ndim,theta_num);
    NumericMatrix theta2(ndim,theta_num);
    NumericMatrix count2(1500,1);
    
    for (int i = 0; i<1500; i++) {
        count2(i,0) = 0;
    }
    
    // set the initial value for mu and theta
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<numbers; j++) {
            mu(i,j) = data(i,j);
            mu2(i,j) = mu(i,j);
        }
    }
    
    for(int i = 0;i<ndim;i++){
        for (int j =0; j<theta_num; j++) {
            int row = j /numbers;
            int col = j %numbers;
            theta(i,j) = mu(i,row)-mu(i,col);
            theta2(i,j) = mu(i,row)-mu(i,col);
        }
    }
    
    double lambda1 = rho;
    NumericVector group(numbers);
    double S1 = 0;
    double S2 = -1;
    double epsilon = 0.001;
    //stopping_criteria(mu,mu1, ndim, numbers, count,epsilon)==0
    while ((S2 < S1 || count < 2 ) && count < 500 && stopping_criteria(mu,mu2, ndim, numbers, count,epsilon)==0) {
        int tempcount = 0;
        
        for(int i = 0;i<ndim;i++){
            for (int j =0; j<numbers; j++) {
                mu2(i,j) = mu(i,j);
            }
        }
        
        while (stopping_criteria2(mu, theta, theta1,u, lambda1, abs_res, rel_res, ndim, numbers, tempcount)==0) {
            if (mumethod ==0) {
            for (int j=0; j<numbers; j++) {
                for (int i = 0; i<ndim; i++) {
                    double temp_sum1=0;
                    double temp_sum2=0;
                    for (int k = (j+1); k < numbers;k++)
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
                    mu(i,j) =  (data(i,j)+lambda1*temp_sum1+lambda1*temp_sum2) / (1 + lambda1*(numbers-1));
                }
            }
            } else {
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
            if (methods ==0 ) { //gTLP
                for (int j=0; j<theta_num; j++) {
                    int row = j / numbers;
                    int col = j % numbers;
                    // only update the theta we need.
                    if (row <col) {
                        if (distance_2(theta2, ndim, j)>= tau) {
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
            } else {
                for (int j=0; j<theta_num; j++) {
                    int row = j / numbers;
                    int col = j % numbers;
                    // only update the theta we need.
                    if (row <col) {
                        for (int i = 0; i<ndim; i++) {
                            double abstheta = theta2(i,j);
                            if (theta2(i,j) <0) {
                                abstheta = -theta2(i,j);
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
            
            // update u
            for (int j=0; j<theta_num; j++) {
                int row = j / numbers;
                int col = j % numbers;
                // only update the theta we need.
                if (row <col) {
                    for(int i = 0;i <ndim;i++)
                    {
                        u(i,j) = u(i,j) + theta(i,j)-(mu(i,row)-mu(i,col));
                    }
                }
            }
            
            // count2(tempcount,0)  = cal_primalRes(mu,theta,ndim,numbers);
            
            tempcount ++;
        }
        
        if (count ==0 ) {
            S1 = S2 = cal_S_ADMM(data,theta, theta2, mu,lambda2,tau,ndim,numbers);
        } else {
            S1 = S2;
            S2 = cal_S_ADMM(data,theta, theta2, mu,lambda2,tau,ndim,numbers);
        }
        
        count2(count,0) =tempcount;
        for(int i = 0;i<ndim;i++){
            for (int j =0; j<theta_num; j++) {
                theta2(i,j) = theta(i,j);
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
    out["count2"] = count2;
    
    out["group"] = group;
    return out;
}
