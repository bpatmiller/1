package QuantSim;

import org.ejml.data.CMatrixRMaj;
import org.ejml.dense.row.CommonOps_CDRM;

public class MeasurementOperation extends SpecialCMatrixRMajLink {
    final int[] bits;

    public MeasurementOperation(int[] bits) {
        this.bits = bits;
    }

    @Override
    public void compute(CMatrixRMaj stateVector) {
        Measure.measure(stateVector, bits);
    }

    @Override
    public CMatrixRMaj getMatrix() {
        return null;
    }
}
