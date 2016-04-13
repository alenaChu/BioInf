import java.util.Scanner;
import java.util.Stack;

/**
 * Created by Chuikina Alena on 16.3.16.
 */
public class Approx {
    public static void main(String[] args) {
        // put your code here
        Scanner scanner = new Scanner(System.in);
        String text1 = scanner.next();
        String text2 = scanner.next();
        int d = new Integer(scanner.next());
        Approx hd = new Approx();

        //   System.out.println(hd.hammingDistance(text2,text1));
        System.out.println(hd.approximatePatternCount(text2,text1,d));
        //    for (Object i: hd.approximatePatternOccurency(text1, text2, d))
        //      System.out.print(i + " ");
    }


    /** count pattern with mismatches in genome
     * @param text
     * @param pattern
     * @param d
     * @return integer = how many times does the pattern occures in text
     */
     public int approximatePatternCount(String text, String pattern, int d){
         int k = pattern.length();
         int count =0;
         for (int i=0; i<= text.length() - k; i++)
             if (hammingDistance(pattern, text.substring(i, i + k)) <= d) count++;
         return count;
     }

    /**  Compute the Hamming distance between two strings
     * @param text1 - two strings of the same size
     * @param text2
     * @return integre distance - number of unequal positions
     */
    public int hammingDistance(String text1,String text2){
        if (text1.length() != text2.length())
                throw new IllegalArgumentException("Strings must be of equal sizes");
        int dist = 0;
        char [] array1 = text1.toCharArray();
        char [] array2 = text2.toCharArray();
        for (int i=0; i< array1.length; i++)
            if (array1[i] != array2[i]) dist++;
        return dist;
    }

    /**  Find all approximate occurrences of a pattern in a string.
     * @param pattern
     * @param text
     * @param d - allowed hamming distance between pattern and text-subtring
     * @return All starting positions where Pattern appears as a substring of Text with at most d mismatches.
     */
    public Iterable approximatePatternOccurency(String pattern, String text, int d){
        Stack<Integer> index= new Stack<>();
        int k = pattern.length();
        for (int i=0; i<= text.length() - k; i++)
            if (hammingDistance(pattern,text.substring(i, i+k)) <= d) index.add(new Integer(i));
        return index;
    }

    /**  Distance between pattern and strings
     * @param pattern
     * @param dna - array
     * @return d(pattern, dna)
     */
    public int d(String pattern, String[] dna){
        int dist = 0;
        for (Object text: dna)
            dist+= d(pattern, (String)text);
        return dist;
    }
    //distance between pattern and one string
    private int d(String pattern, String dnaI){
        int min = Integer.MAX_VALUE;
        int k = pattern.length();
        Approx ap = new Approx();
        for (int i=0; i <= dnaI.length() - k; i++){
            String pattern1 = dnaI.substring(i, i + k);
            int dist = ap.hammingDistance(pattern, pattern1);
            if (dist < min) min = dist;
        }
        return min;
    }
}
