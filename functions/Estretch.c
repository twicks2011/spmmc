double Estretch(double x[][500], int N, double sigma, int stretchStart, int stretchFinish){

  return -sigma*((x[0][stretchFinish]-x[0][stretchStart])*(x[0][stretchFinish]-x[0][stretchStart]) +
		 (x[1][stretchFinish]-x[1][stretchStart])*(x[1][stretchFinish]-x[1][stretchStart]) +
		 (x[2][stretchFinish]-x[2][stretchStart])*(x[2][stretchFinish]-x[2][stretchStart]));
}
