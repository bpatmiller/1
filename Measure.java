package QuantSim;

import org.ejml.data.CMatrixRMaj;
import org.ejml.data.Complex_F32;
import org.ejml.dense.row.CommonOps_CDRM;

import java.lang.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Random;
import java.io.IOException;


class Measure {
    public static void main(String[] args) {
        MathContext mc = new MathContext(100);
        BigDecimal mean = new BigDecimal(0,mc);
        AnuRandom random = new AnuRandom(4096);
        byte[] bytes = random.getBytes();
        for (int j = 0;  j< 1024; j++){
            int val = 0;
            String x = "";
            for (int i = 4 * j; i < 4; i++) {
                x += Integer.toString(bytes[i], 2);
            }
            Random rand = new Random(val);
            mean = mean.add(new BigDecimal(rand.nextDouble(), mc), mc);
        }
        mean  = mean.divide(new BigDecimal(1024, mc), mc);
        System.out.println();
        System.out.println(mean.toString());
    }

    public static int getSeed() {
        //Gets bytes from server, will exit if server inaccessible
        AnuRandom random = new AnuRandom(4);
        String temp = new String(random.getBytes());
        //System.out.println(temp);
        //Gets bytes from server, throws catchable exception if server inaccessible
        try {
            temp = new String(random.getBytesSafe());
        } catch (IOException e) {
            //Handle inaccessible server
        }

        int val = 0;
        for (int i = 0; i < 4; i++) {
            val = val << 8;
            val = val | (random.getBytes()[i] & 0xFF);
        }
        return val;
    }

    public static CMatrixRMaj measure(CMatrixRMaj amplitudes, int[] bits) {
        // get seed
        int seed = getSeed();
        Random rand = new Random(seed);

        for (int index = 0; index < bits.length; index++) {
            // calculate probability that |x> = |1>
            double prob = 0.0;
            int ampsize = amplitudes.numRows;
            int bit = bits[index];

            // statement inspired by a logarithm identity
            if (!(0 <= bit && bit < Math.log(ampsize) / Math.log(2)))
                return null;

            // precalculate some numbers for iteration
            int twobit = ampsize / (int) Math.pow(2, bit);
            int twobit1 = ampsize / (int) Math.pow(2, bit + 1);

            Complex_F32 out = new Complex_F32(0, 0);
            for (int pos = ampsize - 1; pos > 0; pos -= twobit)
                for (int i = 0; i < twobit1; i++) {
                    amplitudes.get(pos - i, 0, out);
                    prob += out.getMagnitude2();
                }
            System.out.println(prob);
            //System.out.print("prob:"+prob+" // ");
            // "measure" |x> with random quantum seed
            int state;
            int x =rand.nextInt(100);
            System.out.println(x);
            if (prob * 100 > x -10) state = 1;
            else state = 0;
            //System.out.println("state:"+state);
            // given the measurement, 0 out entries
            for (int j = (state == 1 ? 0 : twobit1); j < ampsize; j += twobit) {
                for (int k = 0; k < twobit1; ++k)
                    amplitudes.set(k + j, 0, 0, 0);
            }
            // renormalize new state vectors
            float mag = 0;
            for (int i = 0; i < ampsize; i++) {
                amplitudes.get(i, 0, out);
                mag += out.getMagnitude2();
            }
            mag = (float) (1 / Math.sqrt(mag));
            CommonOps_CDRM.scale(mag, 0, amplitudes);
        }
        return amplitudes;
    }
}