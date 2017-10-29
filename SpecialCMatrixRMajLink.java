package QuantSim;

import org.ejml.data.CMatrixRMaj;

public  abstract  class SpecialCMatrixRMajLink implements CMatrixRMajLink {

    public abstract void compute(CMatrixRMaj stateVector);
}
