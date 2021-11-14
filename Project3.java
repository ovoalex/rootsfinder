public class Project3 {

    /**
     * Bisection method.
     */
    private static void bisection(Function func, double a, double b, int maxIterations, double ERROR) {
        System.out.println("BISECTION METHOD:");
        System.out.println("-----------------");

        //find a and b starting point
        double fa = func.func(a);
        double fb = func.func(b);
        if (fa * fb >= 0) {
            System.out.println("Function has same sign at a and b.");
            return;
        }

        double prev_c = 0;
        double prevError = 0;
        int divergeTracker = 0;

        //The maximum iterations stops if the solution is very slowly convergent.
        //Estimate the root by evaluating f(c)
        System.out.printf("%3s %9s %9s %9s %9s %9s %9s %9s%n", "n", "a_n", "b_n", "c_n", "f(a_n)", "f(b_n)", "f(c_n)", "error");
        for (int count = 0; count <= maxIterations; count++) {
        	//Construct the midpoint of the interval c = (a+b)/2.
            double c = (a + b) / 2;
            double fc = func.func(c);
            double error = Math.abs((c - prev_c) / c);

            if (count == 0)
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9s%n", count, a, b, c, fa, fb, fc, "N/A");
            else
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f%n", count, a, b, c, fa, fb, fc, error);
            
            //ERROR desired relative error.
            if (fc == 0.0001 || error < ERROR) {
                System.out.println();
                return;
            }

            // Divergence checker
            if (error > prevError)
                divergeTracker += 1;
            if (divergeTracker > 3) {
                System.out.println("Divergence.");
                return;
            }

            if (fc * fa < 0) {
                b = c;
                fb = fc;
            }
            else {
                a = c;
                fa = fc;
            }
            prev_c = c;
            prevError = error;
        }
    }


    /**
     * False-Position method
     */
    private static void false_position(Function func, double a, double b, int maxIterations, double ERROR) {
        System.out.println("\nFALSE POSITION METHOD:");
        System.out.println("------------------------");
        
        //find a and b starting point
        double fa = func.func(a);
        double fb = func.func(b);
        if (fa *fb >= 0) {
            System.out.println("Function has same sign at a and b.");
            return;
        }

        double prev_c = 0;

        //The maximum iterations stops if the solution is very slowly convergent.
        System.out.printf("%3s %9s %9s %9s %9s %9s %9s %9s%n", "n", "a_n", "b_n", "c_n", "f(a_n)", "f(b_n)", "f(c_n)", "error");
        for (int count = 0; count <= maxIterations; count++) {
        	//Construct the point c = (a * f(b) - b * f(a))/(f(b) - f(a))
        	//Estimate the root by evaluating f(c).
            double c = (a * fb - b * fa) / (fb - fa);
            double fc = func.func(c);
            double error = Math.abs((c - prev_c) / c);
            double prevError = 0;
            int divergeTracker = 0;

            if (count == 0)
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9s%n", count, a, b, c, fa, fb, fc, "N/A");
            else
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f%n", count, a, b, c, fa, fb, fc, error);

            if (fc == 0.00 || error < ERROR) {
                System.out.println();
                return;
            }

            // Divergence checker
            if (error > prevError)
                divergeTracker += 1;
            if (divergeTracker > 3) {
                System.out.println("Divergence.");
                return;
            }

            if (fc * fa < 0) {
                b = c;
                fb = fc;
            }
            else {
                a = c;
                fa = fc;
            }
            prev_c = c;
            prevError = error;
        }
    }


    /**
     * Newton-Raphson method.
     * Construct x_(i+1) = x_i - (f(x_i)) / (f'(x_i))
     * Estimate the root by evaluating f(x_(i+1)).
     */
    public static void newton_raphson(Function func, double x_init, int maxIterations, double ERROR) {
        System.out.println("\nNEWTON-RAPHSON METHOD:");
        System.out.println("------------------------");
        //Function needs to find the root.
        double x_n = x_init;
        double fx_n = func.func(x_n);
        double fx_n_derv = func.func_der(x_n);
        double prevError = 0;
        int divergeTracker = 0;

        System.out.printf("%3s %9s %9s %9s %9s %9s %9s%n", "n", "x_n", "fx_n", "f'(x_n)", "x_n+1", "f(x_n+1)", "error");
        for (int count = 1; count <= maxIterations; count++) {
        	// a The starting point.
            // b The starting point.
        	//The maximum iterations stops if the solution is very slowly convergent.
            double x_n1 = x_n - (fx_n/ fx_n_derv);
            double fx_n1 = func.func(x_n1);
            double error = Math.abs((x_n1 - x_n) / x_n1);

            System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f%n", count, x_n, fx_n, fx_n_derv, x_n1, fx_n1, error);

            if (fx_n1 == 0.00 || error < ERROR) {
                System.out.println();
                return;
            }

            // Divergence checker
            if (error > prevError)
                divergeTracker += 1;
            if (divergeTracker > 3) {
                System.out.println("Divergence.");
                return;
            }

            x_n = x_n1;
            fx_n = fx_n1;
            fx_n_derv = func.func_der(x_n);
            prevError = error;
        }
    }


    /**
     * Secant method.
     * Construct the point  x_(i+1) = x_i - ((x_i - x_(i-1))/(f(x_i) - f(x_(i-1)) * f(x_i).
     * Estimate the root by evaluating f(x_(i+1)).
     */
    public static void secant (Function func, double a, double b, int maxIterations, double ERROR) {
        System.out.println("\nSECANT METHOD:");
        System.out.println("----------------");
        //function needs to find the root.
        double fa = func.func(a);
        double fb = func.func(b);
        double prev_c = 0;
        double prevError = 0;
        int divergeTracker = 0;
       
        System.out.printf("%3s %9s %9s %9s %9s %9s %10s %8s%n", "n", "x_n-1", "x_n", "f(x_n-1)", "f(x_n)", "x_n+1", "f(x_n+1)", "error");
        // a The starting point x_(i-1).
        // b The starting point x_i.
        //The maximum iterations stops if the solution is very slowly convergent.
        for (int count = 0; count <= maxIterations; count++) {
            double c = b - fb * ((b - a) / (fb - fa));
            double fc = func.func(c);
            double error = Math.abs((c - prev_c) / c);

            if (count == 0)
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %10.4f %8s%n", count, a, b, fa, fb, c, fc, "N/A");
            else
                System.out.printf("%3d %9.4f %9.4f %9.4f %9.4f %9.4f %10.4f %8.4f%n", count, a, b, fa, fb, c, fc, error);

            if (fc == 0.00 || error < ERROR) {
                System.out.println();
                return;
            }

            // Divergence checker
            if (error > prevError)
                divergeTracker += 1;
            if (divergeTracker > 3) {
                System.out.println("Divergence.");
                return;
            }

            a = b;
            fa = fb;
            b = c;
            fb = fc;
            prev_c = c;
            prevError = error;
        }
    }


    /**
     * Modified Secant method.
     * Construct the point  x_(i+1) = x_i - f(x_i) * ((δ * x_i)/((f(x_i) + δ * x_i) - f(x_i))).
     * Estimate the root by evaluating f(x_(i+1)).
     */
    public static void modified_secant (Function func, double a, double delta, int maxIterations, double ERROR) {
        System.out.println("\nMODIFIED SECANT METHOD:");
        System.out.println("-------------------------");

        double fa = func.func(a);
        double prev_b = 0;
        double prevError = 0;
        int divergeTracker = 0;
        
       //The maximum iterations stops if the solution is very slowly convergent.
        System.out.printf("%3s %9s %9s %9s %10s %8s%n", "n", "x_n", "x_n+1", "f(x_n)", "f(x_n+1)", "error");
        
        //The function needs to find the root starting point x_i.
        //delta The value for δ is not automatic.
        for (int count = 0; count <= maxIterations; count++) {
            double b = a - fa * ((delta * a) / (func.func(a + delta * a) - fa));
            double fb = func.func(b);
            double error = Math.abs((b - prev_b) / b);

            if (count == 0)
                System.out.printf("%3d %9.4f %9.4f %9.4f %10.4f %8s%n", count, a, b, fa, fb, "N/A");
            else
                System.out.printf("%3d %9.4f %9.4f %9.4f %10.4f %8.4f%n", count, a, b, fa, fb, error);

            
            if (fb == 0.00 || error < ERROR) {
                System.out.println();
                return;
            }

            // Divergence checker
            //Estimate the root by evaluating f(x_(i+1))
            if (error > prevError)
                divergeTracker += 1;
            if (divergeTracker > 3) {
                System.out.println("Divergence.");
                return;
            }

            a = b;
            fa = fb;
            prev_b = b;
            prevError = error;
        }
    }
    

    public static void main (String args []) {
        float ERROR = (float)0.01;
        int maxIterations = 100;
        @SuppressWarnings("unused")
		Project3 findingRoots = new Project3();
        
        // Equation a
        System.out.println("                 **************EQUATION***************");
        System.out.println("                 (a) f(x) = 2x^3 – 11.7x^2 + 17.7x – 5");
        System.out.println("                 *************************************");
        
        Function func1 = new Function1();
        // First root
        System.out.println("First root");
        Project3.bisection(func1, 0, 1, maxIterations, ERROR);
        Project3.false_position(func1, 0, 1, maxIterations, ERROR);
        Project3.newton_raphson(func1, 0.1, maxIterations, ERROR);
        Project3.secant(func1, 0.1, 1, maxIterations, ERROR);
        Project3.modified_secant(func1, 0.1, 0.01, maxIterations, ERROR);
        // Second root
        System.out.println("Second root");
        Project3.bisection(func1, 1, 2, maxIterations, ERROR);
        Project3.false_position(func1, 1, 2, maxIterations, ERROR);
        Project3.newton_raphson(func1, 2.5, maxIterations, ERROR);
        Project3.secant(func1, 1, 2, maxIterations, ERROR);
        Project3.modified_secant(func1, 1.5, 0.01, maxIterations, ERROR);
        // Third root
        System.out.println("Third root");
        Project3.bisection(func1, 3, 4, maxIterations, ERROR);
        Project3.false_position(func1, 3, 4, maxIterations, ERROR);
        Project3.newton_raphson(func1, 4, maxIterations, ERROR);
        Project3.secant(func1, 3, 4, maxIterations, ERROR);
        Project3.modified_secant(func1, 4, 0.01, maxIterations, ERROR);


        // Equation b
        System.out.println("\n");
        System.out.println("                     **************EQUATION*********");
        System.out.println("                     (b) f(x) = x + 10 – xcosh(50/x)");
        System.out.println("                     *******************************");
        Function func2 = new Function2();

        Project3.bisection(func2, 110, 130, maxIterations, ERROR);
        Project3.false_position(func2, 110, 130, maxIterations, ERROR);
        Project3.newton_raphson(func2, 110, maxIterations, ERROR);
        Project3.secant(func2, 110, 130, maxIterations, ERROR);
        Project3.modified_secant(func2, 110, 0.01, maxIterations, ERROR);
    }
}



interface Function {
    double func(double x);
    double func_der(double x);
}
// Equation (a) f(x) = 2x^3 – 11.7x^2 + 17.7x – 5
class Function1 implements Function {
    public double func(double x) {
        return 2*Math.pow(x,3) - 11.7*Math.pow(x,2) + 17.7*x - 5;
    }
    public double func_der(double x) {
        return 6*Math.pow(x,2) - 23.4*x + 17.7;
    }
}

// Equation (b) f(x) = x + 10 – xcosh(50/x)
class Function2 implements Function {
    public double func(double x) {
        return x + 10 - x * Math.cosh(50/x);
    }
    public double func_der(double x) {
        return 1 + ((50 * Math.sinh(50/x)) / x) - Math.cosh(50/x);
    }
}