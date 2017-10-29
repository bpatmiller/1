import java.lang.*;
import java.util.Random;
import java.io.IOException;



class Measure {
  public static void main(String[] args) {
    // test 1
    double[] amp1 = {0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5};
    int[] bit1 = {0, 1, 2};
    measure(amp1, bit1);
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
            val=val<<8;
            val=val|(random.getBytes()[i] & 0xFF);
        }
        return val;
  }

  public static double[] measure(double[] amplitudes, int[] bits) {
    // get seed
    int seed = getSeed();
    Random rand = new Random(seed);

    for (int index = 0; index < bits.length; index++) {
      // calculate probability that |x> = |1>
      double prob = 0.0;
      int ampsize = amplitudes.length;
      int bit = bits[index];

      // statement inspired by a logarithm identity
      if (!(0 <= bit && bit < Math.log(ampsize)/Math.log(2)))
        return null;

      // precalculate some numbers for iteration
      int twobit  = ampsize/(int)Math.pow(2, bit);
      int twobit1 = ampsize/(int)Math.pow(2, bit + 1);

      for (int pos = ampsize-1; pos > 0; pos -= twobit)
        for (int i = 0; i < twobit1 ; i++)
          prob += Math.pow(amplitudes[pos - i], 2);
        //System.out.print("prob:"+prob+" // ");
        // "measure" |x> with random quantum seed
        int state;

        if (prob*100 > rand.nextInt(100)+0.5 ) state=1;
        else state=0;
        //System.out.println("state:"+state);
        // given the measurement, 0 out entries
        for( int j = (state==1 ? 0 : twobit1); j < ampsize; j += twobit ) {
          for (int k = 0; k < twobit1; ++k )
            amplitudes[k+j]=0;
        }
        // renormalize new state vectors
        double mag = 0;
        for (double i : amplitudes) mag+=Math.pow(i,2);
        for (int j = 0; j < ampsize; ++j) amplitudes[j]/=Math.sqrt(mag);
    }
    return amplitudes;
  }
}