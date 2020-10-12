#define EIGEN_USE_MKL_ALL
#include <Rcpp.h>
#include <omp.h>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <cmath>
#include "mrcox_types.h"

// [[Rcpp::depends(RcppEigen)]]
//
void rev_cumsum_assign(const MatrixXd &src, MatrixXd &dest)
{
    const int K = src.cols();
    const int N = src.rows();

    for(int k = 0; k < K; ++k){
        double current = 0;
        for (int i = 0; i < N; ++i){
            current += src(N-1-i, k);
            dest(N-1-i, k) = current;
        }
    }
}

class MCox_aligned
{
    const int N; //Number of observations
    const int K; //Number of responses
    const int p; //Number of predictors
    MapMatd X;
    MapMatd status;
    MapMati rankmin;
    MapMati rankmax;

    std::vector<PermMat> orders; //Permutation matrices such that P_i * X sort X in increasing event time order

    // Store intermediate results
    MatrixXd eta; 
    MatrixXd exp_eta;
    MatrixXd risk_denom;
    MatrixXd outer_accumu;
    MatrixXd residual;

    double get_residual(const MapMatd &v, bool get_val=false){
        //std::cout << "first product\n";
        eta.noalias() = X*v;
        //std::cout << "first product done\n";
        // Get them in the right order
        #pragma omp parallel for
        for(int k = 0; k < K; ++k){
            eta.col(k) = orders[k] * eta.col(k);
        }

        exp_eta.noalias() =  eta.array().exp().matrix();

        #pragma omp parallel for
        for(int k = 0; k < K; ++k){
            // reverse cumsum to get risk_denom
            double current = 0;
            for (int i = 0; i < N; ++i){
                current += exp_eta(N-1-i, k);
                risk_denom(N-1-i, k) = current;
            }

            // adjust ties
            // This won't have aliasing problem
            for(int i = 0; i < N; ++i){
                risk_denom(i, k) = risk_denom(rankmin(i, k), k);
            }
        }
        // get outer accumu
        outer_accumu.noalias() = (status.array()/risk_denom.array()).matrix();
        
        #pragma omp parallel for
        for(int k = 0; k < K; ++k){
            double current = 0;
            for(int i = 0; i < N; ++i){
                current += outer_accumu(i, k);
                outer_accumu(i, k) = current;
            }

            // Adjust ties, this also won't have alias
            for(int i = 0; i < N; ++i){
                outer_accumu(i, k) = outer_accumu(rankmax(i, k), k);
            }
        }
        residual.noalias() = (outer_accumu.array() * exp_eta.array() - status.array()).matrix();

        // Get the order back
        #pragma omp parallel for
        for(int k = 0; k<K; ++k){
            residual.col(k) = orders[k].transpose() * residual.col(k);
        }

        double cox_val = 0;
        if(get_val){
            cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
        }
        return cox_val;
    }

    public:
    MCox_aligned(int N,
                 int K,
                 int p,
                 const double *X,
                 const double *status,
                 const int *rankmin,
                 const int *rankmax,
                 const Rcpp::List order_list) : N(N),
                                                K(K),
                                                p(p),
                                                X(X, N, p),
                                                status(status, N, K),
                                                rankmin(rankmin, N, K),
                                                rankmax(rankmax, N, K),
                                                eta(N, K),
                                                exp_eta(N, K),
                                                risk_denom(N, K),
                                                outer_accumu(N, K),
                                                residual(N, K)
    {
        for (int k = 0; k < K; ++k){
            orders.emplace_back(Rcpp::as<VectorXi>(order_list[k]));
        }
    }



    double get_gradient(const double *vptr, MatrixXd & grad, bool get_val=false){
        MapMatd v(vptr, p, K);
        double cox_val = get_residual(v, get_val);
        //std::cout << "second product\n";
        //grad.noalias() = (residual.transpose() * X).transpose();
        grad.noalias() = X.transpose() * residual;
        //std::cout << "second product done\n";
        return cox_val;
    }

    double get_value_only(const double *vptr){
        MapMatd v(vptr, p, K);
        eta.noalias() = X*v;
        // Get them in the right order
        #pragma omp parallel for
        for(int k = 0; k < K; ++k){
            eta.col(k) = orders[k] * eta.col(k);
        }
        
        exp_eta.noalias() =  eta.array().exp().matrix();

        #pragma omp parallel for
        for(int k = 0; k < K; ++k){
            // reverse cumsum to get risk_denom
            double current = 0;
            for (int i = 0; i < N; ++i){
                current += exp_eta(N-1-i, k);
                risk_denom(N-1-i, k) = current;
            }

            // adjust ties
            // This won't have aliasing problem
            for(int i = 0; i < N; ++i){
                risk_denom(i, k) = risk_denom(rankmin(i, k), k);
            }
        }
        double cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
        return cox_val;
    }

    MatrixXd Rget_residual(const double *vptr){
        MapMatd v(vptr, p, K);
        get_residual(v);
        return residual;
    }
};


void grad_step_l1(Eigen::Map<Eigen::MatrixXd> & B, const MatrixXd & grad, const MatrixXd &v, const double step_size,
                    double lambda_1, const VectorXd & penalty_factor)
{
    B.noalias() = v - step_size*grad;
    B = ((B.cwiseAbs().colwise() - lambda_1*step_size*penalty_factor).array().max(0) * B.array().sign()).matrix();
}

void get_row_norm(const MatrixXd &Bfull, const double step_size, double lambda_2, 
                const VectorXd & penalty_factor, VectorXd & B_row_norm)
{
    B_row_norm.noalias() = Bfull.rowwise().norm().cwiseMax(lambda_2*step_size*penalty_factor);
}

void prox_l2(Eigen::Map<Eigen::MatrixXd> & B, const double step_size, double lambda_2, 
                const VectorXd & penalty_factor, const VectorXd & B_row_norm)
{
    B =  ((B_row_norm.array() - lambda_2*step_size*penalty_factor.array())/(B_row_norm.array())).matrix().asDiagonal() * B;
}

void update_parameters(MatrixXd & B, const MatrixXd & grad, const MatrixXd &v, const double step_size,
                       double lambda_1, double lambda_2, const VectorXd & penalty_factor,
                       VectorXd & B_row_norm)
{
    int K = grad.cols();
    B.noalias() = v - step_size*grad;
    // Apply proximal operator here:
    //Soft-thresholding
    B = ((B.cwiseAbs().colwise() - lambda_1*step_size*penalty_factor).array().max(0) * B.array().sign()).matrix();
    // Group soft-thresholding
    // should be called the pmax of B_row_norm  and lambda_2*step_size
    B_row_norm.noalias() = B.rowwise().norm().cwiseMax(lambda_2*step_size*penalty_factor);
    B =  ((B_row_norm.array() - lambda_2*step_size*penalty_factor.array())/(B_row_norm.array())).matrix().asDiagonal() * B;
}


// [[Rcpp::export]]
Rcpp::List fit_aligned(Rcpp::NumericMatrix X,
                       Rcpp::NumericMatrix status,
                       Rcpp::IntegerMatrix rankmin,
                       Rcpp::IntegerMatrix rankmax,
                       Rcpp::List order_list,
                       Rcpp::NumericMatrix B0,
                       Rcpp::NumericVector lambda_1_all,
                       Rcpp::NumericVector lambda_2_all,
                       VectorXd pfac,
                       double step_size = 1.0,
                       int niter=2000,
                       double linesearch_beta = 1.1,
                       double eps=1e-5 // convergence criteria
                       )
{
    int N = X.rows();
    int p = X.cols();
    int K = status.cols();
    int K0 = B0.cols();
    // B has K >= K0 columns, the extra columns  will be used for group proximal operator
    MapMatd Bmap(&B0(0,0), p, K0);
    MapMatd vmap(&B0(0,0), p, K);

    MCox_aligned prob(N,
                        K,
                        p,
                        &X(0,0),
                        &status(0,0),
                        &rankmin(0,0),
                        &rankmax(0,0),
                        order_list);

    MatrixXd Bfull(Bmap); // has dimension p by K0
    Eigen::Map<Eigen::MatrixXd> B(Bfull.data(), p, K);
    MatrixXd prev_B(vmap);
    MatrixXd v(vmap);
    MatrixXd grad(p,K);
    MatrixXd grad_ls(p,K);
    VectorXd B_row_norm(p);
    const int nlambda = lambda_1_all.size();

    double step_size_intial = step_size;
    double weight_old, weight_new;
    double rhs_ls;
    struct timeval start, end;

    Rcpp::List result(nlambda);
    Rcpp::List residual_result(nlambda);

    for (int lam_ind = 0; lam_ind < nlambda; ++lam_ind){
        gettimeofday(&start, NULL);

        double lambda_1 = lambda_1_all[lam_ind];
        double lambda_2 = lambda_2_all[lam_ind];
        double step_size = step_size_intial;
        weight_old = 1.0;
        bool stop; // stop line search
        double value_change;
        if(lam_ind > 0){
            v.noalias() = B;
        }


        for (int i = 0; i< niter; i++){
            Rcpp::checkUserInterrupt();
            prev_B.noalias() = B;
            double cox_val = prob.get_gradient(v.data(), grad, true);
            while (true){
                grad_step_l1(B, grad, v, step_size, lambda_1, pfac);
                get_row_norm(Bfull, step_size, lambda_2, pfac, B_row_norm);
                prox_l2(B, step_size, lambda_2, pfac, B_row_norm);


                //update_parameters(B, grad, v, step_size, lambda_1, lambda_2, pfac, B_row_norm);

                double cox_val_next = prob.get_value_only(B.data());

                if(!std::isfinite(cox_val_next)){
                    stop = (step_size < 1e-9);
                } else {
                    rhs_ls = cox_val + (grad.array() * (B - v).array()).sum() + (B-v).squaredNorm()/(2*step_size);
                    stop = (cox_val_next <= rhs_ls);
                }


                if (stop){
                    value_change = abs(cox_val_next - cox_val)/fmax(1.0, abs(cox_val));
                    break;
                }
                step_size /= linesearch_beta;
            }


            if(value_change < 5e-7){
                std::cout << "convergence based on value change reached in " << i <<" iterations\n";
                std::cout << "current step size is " << step_size << std::endl;
                gettimeofday(&end, NULL);
                double delta  = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                std::cout <<  "elapsed time is " << delta << " seconds" << std::endl;
                Rcpp::checkUserInterrupt();
                break;
            }

            // Nesterov weight
            weight_new = 0.5*(1+sqrt(1+4*weight_old*weight_old));
            v.noalias() = B + ((weight_old - 1)/weight_new) * (B - prev_B);
            weight_old = weight_new;

            if (i != 0 && i % 100 == 0){
                std::cout << "reached " << i << " iterations\n";
                gettimeofday(&end, NULL);
                double delta  = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                std::cout <<  "elapsed time is " << delta  << " seconds" << std::endl;
                std::cout << "current step size is " << step_size << std::endl;
                Rcpp::checkUserInterrupt();
            }
        }
        result[lam_ind] = Bfull;
        residual_result[lam_ind] = prob.Rget_residual(B.data());
        std::cout << "Solution for the " <<  lam_ind+1 << "th lambda pair is obtained\n";
    }
    return Rcpp::List::create(Rcpp::Named("result") = result,
                              Rcpp::Named("residual") = residual_result);
}


//' @export
// [[Rcpp::export]]
VectorXd compute_dual_norm(MatrixXd grad,
                           double alpha,
                           double tol)
{
    int p = grad.rows();
    VectorXd upperbound((grad.cwiseAbs().rowwise().maxCoeff()).cwiseMin(grad.rowwise().norm()/alpha));
    VectorXd dual_norm(p);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < p; ++i){
        double lower = 0.0;
        double upper = upperbound[i];
        if (upper <= tol){
            dual_norm[i] = 0.0;
        } else {
            int num_iter = (int)ceil(log2(upper/tol));
            for (int j = 0; j < num_iter; ++j){
                double bound = (lower + upper)/2;
                bool less = ((grad.row(i).array().abs() - bound).max(0).matrix().norm()) <= alpha * bound;
                if (less){
                    upper = bound;
                } else {
                    lower = bound;
                }
            }
            dual_norm[i] = (lower + upper)/2;
        }
    }
    return dual_norm;
}

//' @export
// [[Rcpp::export]]
MatrixXd compute_residual(Rcpp::NumericMatrix X,
                       Rcpp::NumericMatrix status,
                       Rcpp::IntegerMatrix rankmin,
                       Rcpp::IntegerMatrix rankmax,
                       Rcpp::List order_list,
                       MatrixXd v)
{
    int N = X.rows();
    int p = X.cols();
    int K = status.cols();

    MCox_aligned prob(N,
                        K,
                        p,
                        &X(0,0),
                        &status(0,0),
                        &rankmin(0,0),
                        &rankmax(0,0),
                        order_list);
    return prob.Rget_residual(v.data());

}