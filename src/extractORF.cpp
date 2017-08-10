#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


// LS: This decleration seems necessary since this function is also 
// used inside extractSeq()
//std::string reverseComplement( std::string );



// [[Rcpp::export]]
SEXP extractSeq(SEXP Gseq, SEXP Left, SEXP Right, SEXP Strand) {
  std::string gseq                = as<std::string>(Gseq);
  std::vector< int > left         = as<std::vector< int > >(Left);
  std::vector< int > right        = as<std::vector< int > >(Right);
  std::vector< int > strand       = as<std::vector< int > >(Strand);
  
  std::vector< std::string > orfs ( left.size() );
  
  for( unsigned i=0; i<orfs.size() ; i++){
    if( left[i] < right[i] ){
      orfs[i] = gseq.substr( left[i]-1, right[i]-left[i]+1 );
    } else {
      orfs[i] = gseq.substr( left[i]-1, gseq.length() );
      orfs[i].append( gseq.substr( 0, right[i] ) );
    }
  }

  return wrap(orfs);
}


// // [[Rcpp::export]]
// std::string reverseComplement( std::string seq ){
//   std::string rcseq = seq;
//   char c;
//   unsigned N = seq.size();
//   
//   for( unsigned i=0; i<N; i++ ){
//     c = (char)seq[i];
//     if( (c == 'A') | (c == 'a') ){
//       rcseq[N-i-1] = 'T';
//     } else if( (c == 'T') | (c == 't') ){
//       rcseq[N-i-1] = 'A';
//     } else if( (c == 'C') | (c == 'c') ){
//       rcseq[N-i-1] = 'G';
//     } else if( (c == 'G') | (c == 'g') ){
//       rcseq[N-i-1] = 'C';
//     } else if( (c == 'R') | (c == 'r') ){
//       rcseq[N-i-1] = 'Y';
//     } else if( (c == 'Y') | (c == 'y') ){
//       rcseq[N-i-1] = 'R';
//     } else {
//       rcseq[N-i-1] = seq[i];
//     }
//   }
//   
//   return rcseq;
// }
