/**
 * Created by Chuikina Alena on 14.3.16.
 */
public class PatternTransformation {

    public static void main(String args[]){
        PatternTransformation pt = new PatternTransformation();
        System.out.println(pt.patternToNumber("TTCGACATGCGCTTAAAATT"));
        System.out.print(pt.numberToPattern(6065, 10));
    }

    /** transform word to number in order to save corresponding to this word info in an array with easy access
     * @param word
     * @return long number
     */
    public long patternToNumber(String word){
        long number = 0;
        int k = word.length();
        for (int i = 0; i < k; i++){
            char c = word.charAt(i);
            if (c == 'A') number = 4*number;
            if (c == 'T') number = 4*number + 3;
            if (c == 'G') number = 4*number + 2;
            if (c == 'C') number = 4*number + 1;
        }
        return number;
    }

    /** decoding of the word
     * @param number
     * @param k - size of encrypted word
     * @return word
     */
    public String numberToPattern(int number, int k){
        StringBuilder word = new StringBuilder();
        for (int i=0; i< k; i++){
            if (number % 4 ==0) word.append('A');
            else if (number % 4 == 1) word.append('C');
            else if (number % 4 == 2) word.append('G');
            else if (number % 4 == 3) word.append('T');
            number = number / 4;
        }
        word.reverse();
        return new String(word);
    }
}
