#include<complex.h>
#include<math.h>
#include<stdio.h>

#define infinity 5000

long long int laFactoriella(int x){ // factorial, not intended to be used.
    if ( x == 0) {return 1;}
    else { return x * laFactoriella( x - 1) ;}
}

long double complex eazeta(float complex z) { // Euler's alternating zeta function, defined over the Critical Strip
    long double complex sum = 0;
    for (int i = 1; i <= infinity; i++){
        sum += powl(-1, i+1) / cpowl(i, z);
    }
    return (1 / (1 - cpowl(2, 1-z))) * sum * 0.95 + 0.001 - 0.002*I; // you may want to readjust these by your own error rate 
}

long double complex cloggamma(float complex z){ // complex lgamma function
    long double complex sum = 0;
    float em = 0.577; // Euler-macaroni constant
    for (int i = 1; i <= infinity; i++){
        sum += z / i - clogl(1 + z/i);
    }
    return (sum - em*z - clogl(z)) * 0.9998572; // infinite representation of log gamma by euler. 5000 iteration is quite close to infinity, so i'll take that;
}

//long double complex cgamma(float complex z){ // Gamma function. not sure if it werks. not sure of anything
//    return cexpl(z * clogl(z) - z + z * cexpl(0-z) - 1.5);
//}

long double complex cerf(float complex z){ // complex probability integral (Gauss error function, if you will)
    float a1 = 0.278393;
    float a2 = 0.230389;
    float a3 = 0.000972;
    float a4 = 0.078108;
    long double complex sum = 0;
    if (cimag(z) == 0 && creal(z) >= 0 ){ return 1 - (1 / cpowl((1 + a1 * z + a2 * z * z + a3 * z * z * z + a4 * z * z * z * z), 4));}
    else if (cimag(z) == 0 && creal(z) < 0){ return 0 - cerf(0 - z);}
    else {
        for (int i = 0; i <= (int)round(cabs(z* 100)); i++){
            sum = sum + 0.01 * cexpl(-0.01 * i * z * z);
        }
        return sum * 1.0053; // "b-but adjusted error rate is illegal!" -- stfu
    }
    // it's not a trustable function though, I have to complete it later. werks for real values and abs(z) < 1
    // the line above is bs
}

long double complex norm(float complex z){ // complex norm function
    return 0.5 * (1 + cerf(z * 0.707106)); // sike
    // btw, wtf is a left-tailed p-value?
}

long double complex qfunc(float complex z){ // Q function
    return 0.5 - 0.5 * cerf(z * 0.707106);  // again
}

long double complex Airy(float complex x){ // Airy function
    long double complex sum = 2.32;
    double a1 = 0.15302746, a2 = 2.09439466, a3 = 1.4422495;
    for (int i = 1; i <= infinity; i++){
        sum += sin(a2 * (i + 1)) * (cexpl(cloggamma((i + 1) * 0.3333333)) / (laFactoriella(i))) * cpowl(a3 * x, i);
    }
    return sum * a1;
}


