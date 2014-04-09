#include "StateSpace.h"

using namespace Rcpp;

class Kalman{
  
  public:
	
	
	//Constructor
    Kalman(NumericMatrix x, NumericVector xZt, NumericMatrix xTt, NumericMatrix xRt, 
          NumericMatrix xHt, NumericMatrix xQt){//system matrices must be set first
      //System matrices
	  Zt_=as<arma::mat>(xZt);
	  Tt_=as<arma::mat>(xTt);
	  Rt=as<arma::mat>(xRt);
	  //Data matrices
      yt=as<arma::mat>(x);
      Ht=as<arma::mat>(xHt);
      Qt=as<arma::mat>(xQt);
      //Dimension scalars
      n=yt.n_cols;
      m=Tt_.n_rows;
      p=Zt_.n_rows;
      r=Qt.n_rows;
      //Filtering matrices
      at.set_size(m,n+1);
      Pt.set_size(m,m*(n+1));
      eps_t.set_size(p,n);
      eta_t.set_size(r,n);
      vt.set_size(p,n);
      Ft.set_size(p,p*n);
      Kt.set_size(m,p*n);
      Lt.set_size(m,m*n);
      //Smoothing matrices
      rt.set_size(m,n+1);
      S_at.set_size(m,n);
      Nt.set_size(m,m*(n+1));
      S_Pt.set_size(m,m*n);
      ut.set_size(p,n);
      Dt.set_size(p,p*n);
      S_eps_t.set_size(p,n);
      S_eta_t.set_size(r,n);
      Svar_eps_t.set_size(p,p*n);
      Svar_eta_t.set_size(r,r*n);
      
      //Time varying system matrices
      if(Zt_.n_cols != m){
            bZt = true;
            Zt.set_size(p, m);  }
      else{ Zt = Zt_;
            bZt = false;  }
      if(Tt_.n_cols != m){
            bTt = true;
            Tt.set_size(m, m);  }
      else{ Tt = Tt_;
            bTt = false;  }

      //Initial values
      at.col(0).fill(0.0);
      Pt.cols(0,(m-1))=arma::diagmat(arma::ones<arma::vec>(m)*1000000);
      rt.col(n).fill(0.0);
      Nt.cols(m*n,m*(n+1)-1).fill(0.0);
      
    }
    
	//--------------------------------------------------
    void kalmanFil(){
      
      LogL=0;
      int tm1,tm2,tp1,tp2;
	  double detFt;
      arma::mat invFt;
	  arma::vec eigVal;
      
        for(int t=0;t<n;t++){
          tm1=t*m;tm2=(t+1)*m-1;
          tp1=t*p;tp2=(t+1)*p-1;
          
          if(bZt)Zt = Zt_.cols(tm1, tm2);
          if(bTt)Tt = Tt_.cols(tm1, tm2);
          
          vt.col(t)=yt.col(t)-Zt*at.col(t);
          Ft.cols(tp1,tp2)=Zt*Pt.cols(tm1,tm2)*Zt.t()+Ht;
		  
		  //Check for positive semidefinite and singularity in Ft, and else set LogL very low and exit loop
		  detFt=det(Ft.cols(tp1,tp2));
		  eigVal = arma::eig_sym(Ft.cols(tp1,tp2));
		  if(eigVal.min()<0 || detFt==0){
			if(t<(m*2)){ //only allow this for the 10 irretations
			//Rcout << "Ft not invertible, determinant of Ft set to 1 and generalized inverse used at t=" << t << std::endl;
			detFt=1;
			invFt=arma::pinv(Ft.cols(tp1,tp2));
			}else{
			//Rcout << "Bad values, exiting kalman filter!" << std::endl;
			LogL=100000000;
			break;
			}
		  }else{		  
          invFt=arma::inv(Ft.cols(tp1,tp2));
		  }
	
          Kt.cols(tp1,tp2)=Tt*Pt.cols(tm1,tm2)*Zt.t()*
                    invFt;
          Lt.cols(tm1,tm2)=Tt-Kt.cols(tp1,tp2)*Zt;
          
          at.col(t+1)=Tt*at.col(t)+Kt.cols(tp1,tp2)*vt.col(t);
             
          Pt.cols(tm1+m,tm2+m)=Tt*Pt.cols(tm1,tm2)*Lt.cols(tm1,tm2).t()+
                              Rt*Qt*Rt.t();
          if(t>=m){
            LogL+=as_scalar(log(detFt)+ 
                    vt.col(t).t()*invFt*vt.col(t));
          }
        }

        LogL=-(n-m)*(double)p/2*log(2*PI)-0.5*LogL;
      }  
    
	//------------------------------------------------
    void kalmanSmooth(){
		arma::mat invFt;
		int tm1,tm2,tp1,tp2,tr1,tr2;

    //rt and Nt is not syncronized with at and Pt. at(1)<=>rt(0) -- at(t)<==>S_at(t-1)
        for(int t=(n-1);t>=0;t--){
          tm1=t*m;tm2=(t+1)*m-1;
          tp1=t*p;tp2=(t+1)*p-1;
          tr1=t*r;tr2=(t+1)*r-1;
          
          if(bZt)Zt = Zt_.cols(tm1, tm2);
          if(bTt)Tt = Tt_.cols(tm1, tm2);
            
          invFt=arma::inv(Ft.cols(tp1,tp2));
          rt.col(t)=Zt.t()*invFt*vt.col(t)+
                  Lt.cols(tm1,tm2).t()*rt.col(t+1);
          
          Nt.cols(tm1,tm2)=Zt.t()*invFt*Zt+
            Lt.cols(tm1,tm2).t()*Nt.cols(tm1+m,tm2+m)*Lt.cols(tm1,tm2);
          
          S_at.col(t)=at.col(t)+Pt.cols(tm1,tm2)*rt.col(t);
          
          S_Pt.cols(tm1,tm2)=Pt.cols(tm1,tm2)-Pt.cols(tm1,tm2)*
                            Nt.cols(tm1,tm2)*Pt.cols(tm1,tm2);
          
          ut.col(t)=invFt*vt.col(t)-Kt.cols(tp1,tp2).t()*rt.col(t+1);
          
          Dt.cols(tp1,tp2)=invFt+Kt.cols(tp1,tp2).t()*Nt.cols(tm1+m,tm2+m)*
                                        Kt.cols(tp1,tp2);
                                        
          S_eps_t.col(t)=Ht*ut.col(t);
          S_eta_t.col(t)=Qt*Rt.t()*rt.col(t);
          
          Svar_eps_t.cols(tp1,tp2)=Ht-Ht*Dt.cols(tp1,tp2)*Ht;
          Svar_eta_t.cols(tr1,tr2)=Qt-Qt*Rt.t()*Nt.cols(tm1+m,tm2+m)*Rt*Qt;
          }
        
    }
	
	//--------------------------------------------------
    void kalmanFil2(){
    //In case of missing observations or forecasting, univariate times series only  
	//In case of forecasting add number of NA to end of observation vector, missing values should also be replaced with NAs
      LogL=0;
      int tm1,tm2,tp1,tp2;
	  double detFt;
      arma::mat invFt;
	  arma::vec eigVal;
      
        for(int t=0;t<n;t++){
          tm1=t*m;tm2=(t+1)*m-1;
          tp1=t*p;tp2=(t+1)*p-1;
          
          if(bZt)Zt = Zt_.cols(tm1, tm2);
          if(bTt)Tt = Tt_.cols(tm1, tm2);
          
          if(!R_IsNA(yt(0,t))) vt.col(t)=yt.col(t)-Zt*at.col(t);
		  else vt.col(t).fill(0.0);
		  
          Ft.cols(tp1,tp2)=Zt*Pt.cols(tm1,tm2)*Zt.t()+Ht;
          		  
		  //Check for positive semidefinite and singularity in Ft, and else set LogL very low and exit loop
		  detFt=det(Ft.cols(tp1,tp2));
		  eigVal = arma::eig_sym(Ft.cols(tp1,tp2));
		  if(eigVal.min()<0 || detFt==0){
			if(t<(m*2)){ //only allow this for the 10 irretations
			//Rcout << "Ft not invertible, determinant of Ft set to 1 and generalized inverse used at t=" << t << std::endl;
			detFt=1;
			invFt=arma::pinv(Ft.cols(tp1,tp2));
			}else{
			//Rcout << "Bad values, exiting kalman filter!" << std::endl;
			LogL=100000000;
			break;
			}
		  }else{		  
          invFt=arma::inv(Ft.cols(tp1,tp2));
		  }
        
          if(!R_IsNA(yt(0,t))) Kt.cols(tp1,tp2)=Tt*Pt.cols(tm1,tm2)*Zt.t()*invFt;
		  else Kt.cols(tp1, tp2).fill(0.0);
		  
          Lt.cols(tm1,tm2)=Tt-Kt.cols(tp1,tp2)*Zt;
          
          at.col(t+1)=Tt*at.col(t)+Kt.cols(tp1,tp2)*vt.col(t);
             
          Pt.cols(tm1+m,tm2+m)=Tt*Pt.cols(tm1,tm2)*Lt.cols(tm1,tm2).t()+
                              Rt*Qt*Rt.t();
          if(t>=m && !R_IsNA(yt(0,t))){
            LogL+=as_scalar(log(det(Ft.cols(tp1,tp2)))+ 
                    vt.col(t).t()*invFt*vt.col(t));
          }
        }

        LogL=-(n-m)*(double)p/2*log(2*PI)-0.5*LogL;
      }  
	  
    //-------------------------------------------------- 
    double optimFunc(NumericMatrix xHt, NumericMatrix xQt){

      Ht=as<arma::mat>(xHt);
      Qt=as<arma::mat>(xQt);
      
      kalmanFil();
      return(LogL);
    }
    //-------------------------------------------------- 
   
    
 
  //---------- System matrices --------------
  arma::mat Zt;
  arma::mat Zt_;
  arma::mat Tt;
  arma::mat Tt_;
  arma::mat Ht;
  arma::mat Qt;
  arma::mat Rt;
  //--------- Filtering matrices ------------
  arma::mat yt;
  arma::mat at;
  arma::mat Pt;
  arma::mat eps_t;
  arma::mat eta_t;
  arma::mat vt;
  arma::mat Ft;
  arma::mat Kt;
  arma::mat Lt;
  //--------- Smoothing matrices ------------
  arma::mat rt;
  arma::mat Nt;
  arma::mat S_at;
  arma::mat S_Pt;
  arma::mat ut;
  arma::mat Dt;
  arma::mat S_eps_t;
  arma::mat S_eta_t;
  arma::mat Svar_eps_t;
  arma::mat Svar_eta_t;
  //-------- Log likelihood -----------------
  double LogL;  
  int n; //Not much point in changing this alone
  //--------- Misc (private)--------------------------
  private:
  bool bZt, bTt;
  int p,m,r;

};


RCPP_MODULE(kalMod){
    using namespace Rcpp ;
	                  
    class_<Kalman>( "StateSpace" )
	
	//Constructor
	.constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>()
	
	//variables............................................
	
	//System matrices
	.field( "Zt", &Kalman::Zt )
	.field( "Zt_", &Kalman::Zt_ )
	.field( "Tt", &Kalman::Tt )
	.field( "Tt_", &Kalman::Tt_ )
	.field( "Rt", &Kalman::Rt )
	.field( "Ht", &Kalman::Ht )
	.field( "Qt", &Kalman::Qt )
	
	//Filtering matrices	
	.field( "yt", &Kalman::yt )
	.field( "at", &Kalman::at )
	.field( "Pt", &Kalman::Pt )
	.field( "Ft", &Kalman::Ft )
	.field( "vt", &Kalman::vt )
	.field( "Kt", &Kalman::Kt )
	.field( "Lt", &Kalman::Lt )
	
	//Smoothing  matrices
	.field( "rt", &Kalman::rt )
	.field( "Nt", &Kalman::Nt )
	.field( "S_at", &Kalman::S_at )
	.field( "S_Pt", &Kalman::S_Pt )
	.field( "S_eps_t", &Kalman::S_eps_t )
	.field( "S_eta_t", &Kalman::S_eta_t )
	
	//Log likelihood
	.field( "LogL", &Kalman::LogL )
	//n number of observation points
	.field( "n", &Kalman::n )
	
	//Functions..................................................
	.method( "kalmanFil", &Kalman::kalmanFil , "Run kalman filter" )
	.method( "kalmanFil2", &Kalman::kalmanFil2 , "Run kalman filter with missing observation/forecasting" )
	.method( "kalmanSmooth", &Kalman::kalmanSmooth , "Run kalman smoother" )
	.method( "optimFunc", &Kalman::optimFunc, "Function called from R optimization tools")
	
	;
	}
	

























































SEXP rcpparma_hello_world(){
	
	arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
	arma::mat m2 = arma::eye<arma::mat>( 3, 3 ) ;
	                     
	List res ;                                  
	res["mat+mat"] = m1 + 3 * ( m1 + m2 );
	return res ;
}

SEXP abso(SEXP x){
NumericVector vX(x);

std::vector<double> y;
y.resize(vX.size());
std::vector<double> z;
z.resize(vX.size());
 
 for(unsigned int i=0; i<y.size(); i++){
			y[i] = vX[i];
	}

	for(unsigned int i=0; i<y.size(); i++){
			z[i] = abs(y[i]);
	}
		
 
 return wrap(z);
}

