package QuantSim;

import org.ejml.data.*;
import org.ejml.dense.row.CommonOps_CDRM;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_CDRM;
import org.ejml.dense.row.mult.MatrixMatrixMult_CDRM;

import java.util.Arrays;

public class Circuit {

    final CMatrixRMajLink[][] circuitOperations;
    final boolean[] specialOpFlags;
    int timeSteps;
    int numQubits;
    int currentTimeStep;
    CMatrixRMaj stateVector;

    public Circuit(int numQubits, int timeSteps) {
        this.numQubits = numQubits;
        this.timeSteps = timeSteps;
        circuitOperations = new CMatrixRMajLink[timeSteps][numQubits];
        specialOpFlags = new boolean[timeSteps];
        reset();
    }


    public static void main(String[] args){


        Circuit teleport = getTeleportationCircuit();
        teleport.reset();
        teleport.getStateVector().set(0, 0, (float)(Math.sqrt(3)), 0);
        teleport.getStateVector().set(4, 0, 0, (float)(Math.sqrt(2)));
        CommonOps_CDRM.scale((float)(1/ Math.sqrt(5)),0,teleport.getStateVector());
        eval(teleport, false);

        Circuit quantumFourierTransform = new Circuit(10, 1);
        psuedoRandomStateVector(quantumFourierTransform.stateVector);
        quantumFourierTransform.setCircuitComponent(0, 0, CMatrixRMajLink.QFT(1024));
        eval(quantumFourierTransform, false);

        Circuit stressTest = getStressTestCircuit(10, 30);
        psuedoRandomStateVector(quantumFourierTransform.stateVector);
        eval(stressTest, false);


    }

    private static Circuit getStressTestCircuit(int numQubits, int numTimeSteps) {
        Circuit c = new Circuit(numQubits, numTimeSteps);
        CMatrixRMajLink.CommonMatrix[] values = CMatrixRMajLink.CommonMatrix.values();
        for (int i = 0; i < numQubits; i ++){
            for (int j = 0; j < numTimeSteps; j++){
                if(j % 2 == 0){
                    c.setCircuitComponent(i, j, values[(int)(Math.random() * values.length)]);
                }else{
                    if(Math.random() > .5){
                        c.setCircuitComponent(0, j, new MeasurementOperation(new int[]{
                                (int)(Math.random() * numQubits),
                                (int)(Math.random() * numQubits),
                                (int)(Math.random() * numQubits),
                                (int)(Math.random() * numQubits),
                                (int)(Math.random() * numQubits),

                        }));
                    }else {
                        c.setCircuitComponent( 0, j, CMatrixRMajLink.ControlledUnitary(
                                numQubits,
                                (int)(Math.random() * numQubits),
                                (int)(Math.random() * numQubits),
                                values[(int)(Math.random() * values.length)].getMatrix()
                        ));
                    }
                }
            }
        }
        return c;
    }

    private static void psuedoRandomStateVector(CMatrixRMaj amplitudes) {
        for (int i = 0; i < amplitudes.numRows; i++){
            amplitudes.set(i,0, (float)Math.random(), (float)Math.random());
        }
        float mag = 0;
        Complex_F32 out = new Complex_F32(0,0);
        for (int i = 0; i < amplitudes.numRows; i++) {
            amplitudes.get(i, 0, out);
            mag += out.getMagnitude2();
        }
        mag = (float) (1 / Math.sqrt(mag));
        CommonOps_CDRM.scale(mag, 0, amplitudes);
    }

    private static void eval(Circuit circuit, boolean printSteps) {
       while(circuit.currentTimeStep != circuit.timeSteps ){
           circuit.advance();
           if (printSteps){
               circuit.getStateVector().print();
           }
       }
       circuit.getStateVector().print();
    }

    public int getTimeSteps() {
        return timeSteps;
    }

    public int getCurrentTimeStep() {
        return currentTimeStep;
    }

    private static Circuit getTeleportationCircuit() {
        Circuit c = new Circuit(3, 7);

        c.setCircuitComponent(0,0, CMatrixRMajLink.CommonMatrix.I);
        c.setCircuitComponent(1,0, CMatrixRMajLink.CommonMatrix.H);
        c.setCircuitComponent(2,0, CMatrixRMajLink.CommonMatrix.I);

        c.setCircuitComponent(0,1, CMatrixRMajLink.ControlledUnitary(3,1,2,
                CMatrixRMajLink.CommonMatrix.X.getMatrix()));

        c.setCircuitComponent(0,2, CMatrixRMajLink.ControlledUnitary(3,0,1,
                CMatrixRMajLink.CommonMatrix.X.getMatrix()));


        c.setCircuitComponent(0,3, CMatrixRMajLink.CommonMatrix.H);
        c.setCircuitComponent(1,3, CMatrixRMajLink.CommonMatrix.I);
        c.setCircuitComponent(2,3, CMatrixRMajLink.CommonMatrix.I);

        c.setCircuitComponent(0,4, new MeasurementOperation(new int[]{0,1}));

        c.setCircuitComponent(0,5, CMatrixRMajLink.ControlledUnitary(3,1,2,
                CMatrixRMajLink.CommonMatrix.X.getMatrix()));

        c.setCircuitComponent(0,6, CMatrixRMajLink.ControlledUnitary(3,0,2,
                CMatrixRMajLink.CommonMatrix.Z.getMatrix()));

        return c;
    }

    public static CMatrixRMaj kron(CMatrixRMajLink[] circuitOperation) {
        //System.out.println(Arrays.toString(circuitOperation));
        CMatrixRMaj kprod = circuitOperation[0].getMatrix();
        for (int i = 1; i < circuitOperation.length; i++){
            if(circuitOperation[i] == null){
                return kprod;
            }
            CMatrixRMaj nextFactor = circuitOperation[i].getMatrix();
            CMatrixRMaj nextProd = new CMatrixRMaj(kprod.numRows * nextFactor.numRows,
                    kprod.numRows * nextFactor.numRows);
            kron(kprod, nextFactor, nextProd);
            kprod = nextProd;

        }
        return kprod;
    }

    public static void kron(CMatrixRMaj A , CMatrixRMaj B , CMatrixRMaj C )
    {
        int numColsC = A.numCols*B.numCols;
        int numRowsC = A.numRows*B.numRows;

        if( C.numCols != numColsC || C.numRows != numRowsC) {
            throw new IllegalArgumentException("C does not have the expected dimensions");
        }

        // TODO see comment below
        // this will work well for small matrices
        // but an alternative version should be made for large matrices
        Complex_F32 a, b;
        a= new Complex_F32(0,0);
        b = new Complex_F32(0,0);
        for( int i = 0; i < A.numRows; i++ ) {
            for( int j = 0; j < A.numCols; j++ ) {
                A.get(i,j, a);
                for( int rowB = 0; rowB < B.numRows; rowB++ ) {
                    for( int colB = 0; colB < B.numCols; colB++ ) {
                        B.get(rowB, colB, b);
                        Complex_F32 val = a.times(b);
                        C.set(i*B.numRows+rowB,j*B.numCols+colB, val.real, val.imaginary);
                    }
                }
            }
        }
    }

    public CMatrixRMaj getStateVector() {
        return stateVector;
    }

    public void advance(){
        if (currentTimeStep < timeSteps){

            if(specialOpFlags[currentTimeStep]){
                SpecialCMatrixRMajLink specialOp = (SpecialCMatrixRMajLink)circuitOperations[currentTimeStep][0];
                specialOp.compute(stateVector);
                currentTimeStep++;
                return;
            }

            //Recursively compute tensor product of matrices
            //In this row of  circuit operations.
            CMatrixRMaj tensorProduct =kron(circuitOperations[currentTimeStep]);


            //Also probably use a M for measurement

            CMatrixRMaj newStateVector = stateVector.createLike();
            MatrixMatrixMult_CDRM.mult_reorder(tensorProduct, stateVector, newStateVector);
            stateVector = newStateVector;
            currentTimeStep++;
        }
    }

    //Multi bit components should belong to the first bit, and put null in the spaces below.
    public void setCircuitComponent(int bit, int timeStep, CMatrixRMajLink component){
        if(!specialOpFlags[timeStep]){
            circuitOperations[timeStep][bit] = component;
            if(component instanceof  SpecialCMatrixRMajLink){
                specialOpFlags[timeStep] = true;
            }}
    }

    public void reset(){
        float[] stateVectorData = new float[2 * (int)Math.pow(2, numQubits)];
        stateVectorData[0] = 1;
        stateVector = new CMatrixRMaj(stateVectorData.length/2, 1, true, stateVectorData);
        currentTimeStep = 0;
    }
}
