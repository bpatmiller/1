package QuantSim;

import org.ejml.data.CMatrixRMaj;
import org.ejml.data.Complex_F32;
import org.ejml.dense.row.CommonOps_CDRM;

import static QuantSim.Circuit.kron;

public interface CMatrixRMajLink {

     float rt2  =  (float)(Math.sqrt(2));
     float hf = 1/(float)2;
    public enum CommonMatrix implements CMatrixRMajLink{
        I(new CMatrixRMaj(2,2, true,  1,0, 0, 0, 0,0,1,0)),
        X(new CMatrixRMaj(2,2, true,  0,0, 1, 0, 1,0,0,0)),
        H(new CMatrixRMaj(2,2, true,  1/rt2,0, 1/rt2, 0, 1/rt2,0,-1/rt2,0)),
        Y(new CMatrixRMaj(2,2, true,  0,0, 0, -1, 0,1,0,0)),
        Z(new CMatrixRMaj(2,2, true,  1,0, 0, 0, 0,0,-1,0)),
        rtNOT(new CMatrixRMaj(2,2, true, hf, hf , hf,  -hf, hf, -hf, hf, hf));

        private final CMatrixRMaj myMatrix;
        CommonMatrix(CMatrixRMaj myMatrix){
             this.myMatrix = myMatrix;
        }

        @Override
        public CMatrixRMaj getMatrix() {
            return myMatrix;
        }
    }

    public CMatrixRMaj getMatrix();

    public static CMatrixRMajLinkWrapper rotationMatrix(float theta){
        return new CMatrixRMajLinkWrapper(
                new CMatrixRMaj(2,2, true,  1,0, 0, 0, 0,0,
                        (float)Math.cos(theta),(float)Math.sin(theta))
        );
    }

    //n  = 2^#qbits
    public static CMatrixRMajLinkWrapper groverDiffusion(int n){
        float fill = (float)(2/Math.pow(2, n));
        float[] data = new float[2 * n  * n];
        for (int i =0; i < data.length; i+=2){
            data[i] = fill;
        }
        CMatrixRMaj grover = new CMatrixRMaj(n,n, true, data);
        for (int i = 0; i < n; i++){
            grover.set(i,i, fill - 1, 0);
        }
        return new CMatrixRMajLinkWrapper(grover);
    }

    //n = 2^#qbits
    public static CMatrixRMajLinkWrapper QFT(int n){
        Complex_F32 omega0 = new Complex_F32(
                (float)Math.cos(2 * Math.PI/n),
                (float)Math.sin(2 * Math.PI/n)
        );
        //System.out.println(omega0);
        Complex_F32 omegaMult = new Complex_F32(1.0f, 0);
        float[] data = new float[2 * n  * n];
        for (int i = 0; i < n; i++){
            Complex_F32 j_0 = new Complex_F32(1.0f,0);
            for (int j = 2 * n * i; j < 2 * n * (i + 1); j+=2){
                data[j] = j_0.real;
                data[j + 1] = j_0.imaginary;
                j_0 = j_0.times(omegaMult);
            }
            omegaMult = omegaMult.times(omega0);
        }
        CMatrixRMaj qftMatrix = new CMatrixRMaj(n,n, true, data);
        CommonOps_CDRM.scale(1/(float)Math.sqrt(n), 0, qftMatrix);
        return new CMatrixRMajLinkWrapper(qftMatrix);

    }

    //n= 2^#num qbits
    public static CMatrixRMajLinkWrapper ControlledUnitary(int numBits, int control, int target, CMatrixRMaj U){
        float zero = 0;
        float one = 1;
        CMatrixRMajLink[] ops = new CMatrixRMajLink[numBits];

        for (int i = 0; i < numBits; i++)
            ops[i] = CommonMatrix.I;

        float[] control0 = {one, zero, zero, zero, zero, zero, zero, zero};
        ops[control] = new CMatrixRMajLinkWrapper(
                        new CMatrixRMaj(2,2, true,
                        (control0)));
        CMatrixRMaj summand1 = kron(ops);
        float[] control1 = {zero, zero, zero, zero,zero, zero, one, zero};
        ops[control] = new CMatrixRMajLinkWrapper(
                new CMatrixRMaj(2,2, true,
                        (control1)));
        ops[target] = new CMatrixRMajLinkWrapper(U);

        CMatrixRMaj summand2 = kron(ops);
        CMatrixRMaj result = new CMatrixRMaj(summand1.numRows, summand1.numCols);
        CommonOps_CDRM.add(summand1, summand2, result);
        return new CMatrixRMajLinkWrapper(result);
    }


    class CMatrixRMajLinkWrapper implements CMatrixRMajLink{
        private final CMatrixRMaj myMatrix;
        CMatrixRMajLinkWrapper(CMatrixRMaj myMatrix){
            this.myMatrix = myMatrix;
        }

        @Override
        public CMatrixRMaj getMatrix() {
            return myMatrix;
        }
    }
}
