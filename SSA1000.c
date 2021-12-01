// gcc -o ssa_lvpp ssa_lvpp.c
// gcc -O3 -funroll-loops -ftree-vectorize -fast -floop-optimize -fstrength-reduce -o ssa_lvpp ssa_lvpp.c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define MAX_DATA 1000
int main(int argc, char** argv)
{
    // Extract parameters from the command line

    // lambda_e, lambda_f, beta_e, beta_e1, beta_f, beta_f1, d, a, p, 0.0, c, b, g, h, mu


    double    tmax = strtod(argv[1], (char **) NULL);
    double    Xe_start = strtod(argv[2], (char **) NULL);
    double    Ye_start = strtod(argv[3], (char **) NULL);
    double    Xf_start = strtod(argv[4], (char **) NULL);
    double    Yf_start = strtod(argv[5], (char **) NULL);
    double    Ze_start = strtod(argv[6], (char **) NULL);
    double    Zf_start = strtod(argv[7], (char **) NULL);
    double    Y1e_start = strtod(argv[8], (char **) NULL);
    double    Y1f_start = strtod(argv[9], (char **) NULL);

    double    lambdae = strtod(argv[10], (char **) NULL);
    double    lambdaf = strtod(argv[11], (char **) NULL);
    double    beta = strtod(argv[12], (char **) NULL);
    double    beta1 = strtod(argv[13], (char **) NULL);
    double    d = strtod(argv[14], (char **) NULL);
    double    a = strtod(argv[15], (char **) NULL);
    double    p = strtod(argv[16], (char **) NULL);
    double    eta = strtod(argv[17], (char **) NULL);
    double    c = strtod(argv[18], (char **) NULL);
    double    b = strtod(argv[19], (char **) NULL);
    double    g = strtod(argv[20], (char **) NULL);
    double    h = strtod(argv[21], (char **) NULL);
    double    mu = strtod(argv[22], (char **) NULL);

      printf("Xe = %f\n", Xe_start);
      printf("Ye = %f\n", Ye_start);
      printf("Xf = %f\n", Xf_start);
      printf("Yf = %f\n", Yf_start);
      printf("Ze = %f\n", Ze_start);
      printf("Zf = %f\n", Zf_start);

      printf("lambda_e = %f\n", lambdae);
      printf("lambda_f = %f\n", lambdaf);
      printf("beta = %f\n", beta);
      printf("beta1 = %f\n", beta1);
      printf("d = %f\n", d);
      printf("a = %f\n", a);
      printf("p = %f\n", p);
      printf("eta = %f\n", eta);
      printf("c = %f\n", c);
      printf("b = %f\n", b);
      printf("g = %f\n", g);
      printf("h = %f\n", h);
      printf("mu = %f\n", mu);



  //  if (tmax>(MAX_DATA-1)) {printf("tmax=%d\n",MAX_DATA-1); exit(-1);}



    double t = 0;
    double tn = 0;
    int row = 0;
    double data[5];



    FILE *fp = NULL;

    fp = fopen(argv[23] ,"w");


    int run_number = 1;


    double    Xe = Xe_start;
    double    Ye = Ye_start;
    double    Xf = Xf_start;
    double    Yf = Yf_start;
    double    Ze = Ze_start;
    double    Zf = Zf_start;
    double    Y1e = Y1e_start;
    double    Y1f = Y1f_start;

    srand((unsigned)time(NULL));

    // Start the simulation
    for (run_number = 1; run_number < 1001; run_number++){

      //printf("flag1 = %i\n", flag1);
      //printf("flag2 = %i\n", flag2);


      int flag1 = 0;
      int flag2 = 0;

      double first = 0.0;
      double final = 0.0;
      int winner = 0;


      data[0] = run_number;

      t = 0;

      Xe = Xe_start;
      Ye = Ye_start;
      Xf = Xf_start;
      Yf = Yf_start;
      Ze = Ze_start;
      Zf = Zf_start;
      Y1e = Y1e_start;
      Y1f = Y1f_start;


      do{

        // Determine the total rate of each event
        double dXe1 = lambdae;                  // Uninfected EF CD4 created
        double dXe2 = d*Xe;                     // EF CD4 dies naturally
        double dXe3 = beta*Xe*Ye;               // EF CD4 infected by wild type
        double dXe4 = beta1*Xe*Y1e;             // EF CD4 infected by mutant
        double dXe = dXe1 + dXe2 + dXe3 + dXe4;

        double dYe1 = a*Ye;                     // Infected EF CD4 dies
        double dYe2 = p*Ye*Ze;                  // Infected EF CD4 killed by CD8
        double dYe3 = eta*Yf;                   // Infected F CD4 migrates to EF
        double dYe4 = eta*Ye;                   // Infected EF CD4 migrates to F
        double dYe = dYe1 + dYe2 + dYe3 + dYe4;

        double dXf1 = lambdaf;
        double dXf2 = d*Xf;
        double dXf3 = beta*Xf*Yf;
        double dXf4 = beta1*Xf*Y1f;
        double dXf = dXf1 + dXf2 + dXf3 + dXf4;

        double dYf1 = a*Yf;
        double dYf2 = p*Yf*Zf;
        double dYf = dYf1 + dYf2;

        double dZe1 = c*Ye*Ze;
        double dZe2 = b*Ze;
        double dZe3 = g*Ze;
        double dZe4 = h*Zf;
        double dZe = dZe1 + dZe2 + dZe3 + dZe4;

        double dZf = b*Zf;

        double dY1e1 = a*Y1e;
        double dY1e2 = eta*Y1f;
        double dY1e3 = eta*Y1e;
        double dY1e = dY1e1 + dY1e2 + dY1e3;

        double dY1f = a*Y1f;

        // Determine the total event rate
        double toto = dXe + dYe + dXf + dYf + dZe + dZf + dY1e + dY1f;

        // Determine the time step at which the next event occurs + update the time
        double r = (double)rand() / (double)RAND_MAX;

        t += -1/toto * log10(r);

        // Determine the next event to occur + update the populations
        r = (double)rand() / (double)RAND_MAX;
        double r2 = (double)rand() / (double)RAND_MAX;
        double prob = r * toto;


        if (prob <= dXe1) { Xe++;} // Xe ---> Xe + Xe

        else if (prob <= dXe1 + dXe2) {Xe--;} // Xe + Xe ---> Xe

        else if (prob <= dXe1 + dXe2 + dXe3) {
          if (r2 <= mu) {Xe--; Y1e++;} else {Xe--; Ye++;};}   // Xe + Ye ---> Ye + Ye OR Y1e + Ye

        else if (prob <= dXe) {Xe--; Y1e++;}   // Xe + Y1e ---> Y1e + Y1e


        else if (prob <= dXe + dYe1) {Ye--;} // Ye + Ye ---> Ye

        else if (prob <= dXe + dYe1 + dYe2) {Ye--;} // Ye + Ze ---> Ze

        else if (prob <= dXe + dYe1 + dYe2 + dYe3) {Ye++; Yf--;} // Yf ---> Ye

        else if (prob <= dXe + dYe) {Ye--; Yf++;} // Ye ---> Yf


        else if (prob <= dXe + dYe + dXf1) {Xf++;} // Xf ---> Xf + Xf

        else if (prob <= dXe + dYe + dXf1 + dXf2) {Xf--;} // Xf + Xf ---> Xf

        else if (prob <= dXe + dYe + dXf1 + dXf2 + dXf3) {
          if (r2 <= mu) {Xf--; Y1f++;} else {Xf--; Yf++;};} // Xf + Yf ---> Yf + Yf OR Y1f + Yf

        else if (prob <= dXe + dYe + dXf) {Xf--; Y1f++;}    // Xf + Y1f ---> Y1f + Y1f


        else if (prob <= dXe + dYe + dXf + dYf1) {Yf--;}

        else if (prob <= dXe + dYe + dXf + dYf) {Yf--;}


        else if (prob <= dXe + dYe + dXf + dYf + dZe1) {Ze++;}

        else if (prob <= dXe + dYe + dXf + dYf + dZe1 + dZe2) {Ze--;}

        else if (prob <= dXe + dYe + dXf + dYf + dZe1 + dZe2 + dZe3) {Ze--; Zf++;}

        else if (prob <= dXe + dYe + dXf + dYf + dZe) {Ze++; Zf--;}


        else if (prob <= dXe + dYe + dXf + dYf + dZe + dZf) {Zf--;}


        else if (prob <= dXe + dYe + dXf + dYf + dZe + dZf + dY1e1) {Y1e--;}

        else if (prob <= dXe + dYe + dXf + dYf + dZe + dZf + dY1e1 + dY1e2) {Y1e++; Y1f--;}

        else if (prob <= dXe + dYe + dXf + dYf + dZe + dZf + dY1e) {Y1e--; Y1f++;}


        else Y1f--;

        if ((flag1 == 0) && (Y1e > 0)){
          flag1 = 1;
          first = t;
          winner = 1;

          //rintf("%s\n", "Y1e");
        }

        if ((flag1 == 0) && (Y1f > 0)){
          flag1 = 1;
          first = t;
          winner = 2;

          //printf("%s\n", "Y1f");
        }

        if ((flag1 == 1) && (Y1f == 0) && (Y1e == 0)){
          flag1 = 0;
          first = 0.0;
          winner = 0;

          //printf("%s\n", "Y1f");
        }


        if ((flag2 == 0) && (((Y1e + Y1f)/(Y1e + Y1f + Ye + Yf)) >= 0.9)){
          flag2 = 1;
          final = t;

        }


    } while (t<=tmax && (Xe + Ye + Xf + Yf + Ze + Zf + Y1e + Y1f > 0) && flag2 == 0);
    double dTime = final - first;
    data[1] = first;
    data[2] = winner;
    data[3] = final;
    data[4] = dTime;

    if(dTime > 0) {fprintf(fp, "%f \t %f \t %f \t %f \t %f \n", data[0], data[1], data[2], data[3], data[4]);}

    flag1 = 0;
    flag2 = 0;



  }

    fclose(fp);

    exit(EXIT_SUCCESS);
}
