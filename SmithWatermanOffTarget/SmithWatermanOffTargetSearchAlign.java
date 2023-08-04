package SmithWatermanOffTarget;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import java.nio.file.Path;
import java.nio.file.Files;
import java.io.File;
import java.io.InputStream;
import java.io.IOException;

import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;

import java.time.Instant;
import java.time.Duration;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign;

public class SmithWatermanOffTargetSearchAlign {
    private static final String CHROMOSOME_STR = "Chromosome";
    private static final String END_POS_STR = "EndPosition";
    private static final String SEQ_STR = "SiteSeqPlusMaxEditsBefore";
    private static final String EDIT_NUM_STR = "#Edit";
    private static final String STRAND_STR = "Strand";
    private static final String ALINGED_TARGET_STR = "AlignedTarget";
    private static final String ALINGED_TEXT_STR = "AlignedText";
    private static final String MISMATCHES_NUM_STR = "#Mismatches";
    private static final String BUGLES_NUM_STR = "#Bulges";
    private static int NUM_THREADS = 8;
    private static int SITE_WINDOW_SIZE = 50;
    private static int MAX_EDITS = 7;
    private static int MAX_MISMATCHES_WITHOUT_BULGES = MAX_EDITS;
    private static int MAX_MISMATCHES_WITH_BULGES = 6;
    private static int MAX_BULGES = 1;
    private static int EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    private static boolean ALLOW_PAM_EDITS = false;


    public static void setNumThreads(int value) {
        NUM_THREADS = value;
    }

    public static void setSiteWindowSize(int value) {
        SITE_WINDOW_SIZE = value;
    }

    public static void setMaxEdits(int value) {
        MAX_EDITS = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setMaxMismatchesWithoutBulges(int value) {
        MAX_MISMATCHES_WITHOUT_BULGES = value;
    }

    public static void setMaxMismatchesWithBulges(int value) {
        MAX_MISMATCHES_WITH_BULGES = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setMaxBulges(int value) {
        MAX_BULGES = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setAllowPamEdits(boolean value) {
        ALLOW_PAM_EDITS = value;
    }

    public static Alignment find_alignment(int[][] M, String pam, String strand,
        Boolean allowNsInText, int trueDistance, String target, String text, int target_i, int text_i,
        StringBuilder[] targetTextAlign, int mismatches, int bulges) {
        int targetLen = target.length();

        if (bulges > MAX_BULGES) {
            return null;
        }

        if (bulges > 0 && mismatches > MAX_MISMATCHES_WITH_BULGES) {
            return null;
        }
        if (mismatches > MAX_MISMATCHES_WITHOUT_BULGES) {
            return null;
        }
        if (target_i == 0) {
                return new Alignment(mismatches, bulges, targetTextAlign[0], targetTextAlign[1]);
        }
        if (text_i <= 0) {
            return null;
        }
        int effective_max_edit = bulges == 0?  MAX_EDITS : EFFECTIVE_MAX_EDIT_WITH_BULGE;
        if (M[target_i][text_i] > (effective_max_edit - mismatches - bulges)) {
            return null;
        }

        // mismatch
        StringBuilder[] targetTextAlignMissOrMatch = new StringBuilder[2];
        targetTextAlignMissOrMatch[0] = new StringBuilder(targetTextAlign[0].toString());
        targetTextAlignMissOrMatch[1] = new StringBuilder(targetTextAlign[1].toString());
        targetTextAlignMissOrMatch[0].append(target.charAt(target_i - 1));
        targetTextAlignMissOrMatch[1].append(text.charAt(text_i - 1));
        int mismatchSocre;
        if (allowNsInText) {
            mismatchSocre = (target.charAt(target_i - 1) != text.charAt(text_i - 1) &&
                                 target.charAt(target_i - 1) != 'N' && text.charAt(text_i - 1) != 'N') ? 1 : 0;
        }
        else {
            mismatchSocre = (target.charAt(target_i - 1) != text.charAt(text_i - 1) &&
                                 target.charAt(target_i - 1) != 'N') ? 1 : 0;
        }
        if (!ALLOW_PAM_EDITS) {
            if  (((strand == "+") && (target_i <= targetLen &&  targetLen - pam.length() < target_i)) ||
                 ((strand == "-") && (target_i >= 1 &&  target_i <= pam.length()))) {
                    // TODO: shorten this
                    if (strand == "-" && target_i == pam.length()){
                        // This is a special case as we might have "-" before the PAM
                        Alignment targetTextMissOrMatchAlignment;
                        if (mismatchSocre != 0) {
                            targetTextMissOrMatchAlignment = null;
                        }
                        else {
                            targetTextMissOrMatchAlignment = find_alignment(M, pam, strand,
                            allowNsInText, trueDistance, target, text, target_i - 1, text_i - 1,
                            targetTextAlignMissOrMatch, mismatches + mismatchSocre, bulges);
                        }
                        if (targetTextMissOrMatchAlignment != null &&
                            ((targetTextMissOrMatchAlignment.getBulges() + targetTextMissOrMatchAlignment.getMismatches()) == trueDistance)) {
                                return targetTextMissOrMatchAlignment;
                            }
                        StringBuilder[] targetTextAlignTextBugle = new StringBuilder[2];
                        targetTextAlignTextBugle[0] = new StringBuilder(targetTextAlign[0].toString());
                        targetTextAlignTextBugle[1] = new StringBuilder(targetTextAlign[1].toString());
                        targetTextAlignTextBugle[0].append('-');
                        targetTextAlignTextBugle[1].append(text.charAt(text_i - 1));
                        Alignment targetTextnTextBugleAlignment = find_alignment(M, pam, strand,
                            allowNsInText, trueDistance, target, text, target_i, text_i - 1,
                            targetTextAlignTextBugle, mismatches, bulges + 1);
                        if (targetTextnTextBugleAlignment != null &&
                            ((targetTextnTextBugleAlignment.getBulges() + targetTextnTextBugleAlignment.getMismatches()) == trueDistance)) {
                            return targetTextnTextBugleAlignment;
                        }
                        if (targetTextMissOrMatchAlignment == null){
                            return targetTextnTextBugleAlignment;
                        }
                        if (targetTextnTextBugleAlignment == null){
                            return targetTextMissOrMatchAlignment;
                        }
                        if ((targetTextMissOrMatchAlignment.getBulges() + targetTextMissOrMatchAlignment.getMismatches()) <= (targetTextnTextBugleAlignment.getBulges() + targetTextnTextBugleAlignment.getMismatches())){
                            return targetTextMissOrMatchAlignment;
                        }
                        return targetTextnTextBugleAlignment;
                    }
                    
                    if (mismatchSocre != 0) {
                        return null;
                    }
                    return find_alignment(M, pam, strand,
                        allowNsInText, trueDistance, target, text, target_i - 1, text_i - 1,
                        targetTextAlignMissOrMatch, mismatches + mismatchSocre, bulges);
                 }
        }
        Alignment targetTextMissOrMatchAlignment = find_alignment(M, pam, strand,
            allowNsInText, trueDistance, target, text, target_i - 1, text_i - 1, targetTextAlignMissOrMatch, mismatches + mismatchSocre, bulges);
        if (targetTextMissOrMatchAlignment != null &&
            ((targetTextMissOrMatchAlignment.getBulges() + targetTextMissOrMatchAlignment.getMismatches()) == trueDistance)) {
            return targetTextMissOrMatchAlignment;
        }

        // target bulge
        StringBuilder[] targetTextAlignTargetBugle = new StringBuilder[2];
        targetTextAlignTargetBugle[0] = new StringBuilder(targetTextAlign[0].toString());
        targetTextAlignTargetBugle[1] = new StringBuilder(targetTextAlign[1].toString());
        targetTextAlignTargetBugle[0].append(target.charAt(target_i - 1));
        targetTextAlignTargetBugle[1].append('-');
        Alignment targetTextTargetBugleAlignment = find_alignment(M, pam, strand,
            allowNsInText, trueDistance, target, text, target_i - 1, text_i, targetTextAlignTargetBugle, mismatches, bulges + 1);
        if (targetTextTargetBugleAlignment != null &&
            ((targetTextTargetBugleAlignment.getBulges() + targetTextTargetBugleAlignment.getMismatches()) == trueDistance)) {
            return targetTextTargetBugleAlignment;
        }

        // text bulge
        Alignment targetTextnTextBugleAlignment = null;
        if (target_i != targetLen) {
            // we do to have text bulge in the begining/start of the alignment
            // note that it impossible to put text bulge when target_i == 0
            StringBuilder[] targetTextAlignTextBugle = new StringBuilder[2];
            targetTextAlignTextBugle[0] = new StringBuilder(targetTextAlign[0].toString());
            targetTextAlignTextBugle[1] = new StringBuilder(targetTextAlign[1].toString());
            targetTextAlignTextBugle[0].append('-');
            targetTextAlignTextBugle[1].append(text.charAt(text_i - 1));
            targetTextnTextBugleAlignment = find_alignment(M, pam, strand,
                allowNsInText, trueDistance, target, text, target_i, text_i - 1, targetTextAlignTextBugle, mismatches, bulges + 1);
        }
        if (targetTextnTextBugleAlignment != null &&
            ((targetTextnTextBugleAlignment.getBulges() + targetTextnTextBugleAlignment.getMismatches()) == trueDistance)) {
            return targetTextnTextBugleAlignment;
        }
        
        if (targetTextMissOrMatchAlignment != null || targetTextTargetBugleAlignment != null || targetTextnTextBugleAlignment != null) {
            int missOrMatchEdit = targetTextMissOrMatchAlignment != null ?
                targetTextMissOrMatchAlignment.getBulges() + targetTextMissOrMatchAlignment.getMismatches() : MAX_EDITS;
            int targetBugleEdit = targetTextTargetBugleAlignment != null ?
                targetTextTargetBugleAlignment.getBulges() + targetTextTargetBugleAlignment.getMismatches() : MAX_EDITS;
            int textBugleEdit = targetTextnTextBugleAlignment != null ?
            targetTextnTextBugleAlignment.getBulges() + targetTextnTextBugleAlignment.getMismatches() : MAX_EDITS;
            int minEdit = Math.min(missOrMatchEdit, Math.min(targetBugleEdit, textBugleEdit));
            if (targetTextMissOrMatchAlignment!= null && missOrMatchEdit == minEdit) {
                return targetTextMissOrMatchAlignment;
            }
            if (targetTextTargetBugleAlignment!= null && targetBugleEdit == minEdit) {
                return targetTextTargetBugleAlignment;
            } 
            return targetTextnTextBugleAlignment;
        }
        return null;
    }
    
    private static void smithWatermanLastCellLogic(
        int[][] M, String target,String text, int row, int col, Boolean allowNsInText) {
        int matchCell = M[row - 1][col - 1];
        int targetBulgeCell = M[row - 1][col];
        int textBulgeCell = M[row][col - 1];
        if (allowNsInText) {
            M[row][col] = Math.min(Math.min(
                (target.charAt(row - 1) != text.charAt(col - 1) &&
                 target.charAt(row - 1) != 'N' &&
                 text.charAt(col - 1) != 'N') ? matchCell + 1 : matchCell,
                 targetBulgeCell + 1), textBulgeCell + 1
                );
        }
        else {
            M[row][col] = Math.min(Math.min(
                (target.charAt(row - 1) != text.charAt(col - 1) &&
                 target.charAt(row - 1) != 'N') ? matchCell + 1 : matchCell,
                 targetBulgeCell + 1), textBulgeCell + 1
                );
        }
    }

    private static void smithWatermanRowsFill(int[][] M, int m, int n,  String target, String text, Boolean allowNsInText) {
        // fill the rows except to the last row
        for (int row = 1; row <= m - 1; row++) {
            for (int col = 1; col <= n; col++) {
                smithWatermanLastCellLogic(M, target,text, row, col, allowNsInText);
            }
        }
    }

    private static void smithWatermanProcessAlignment(
        Alignment targetTextAlignment, int targetTextAlignmentPos, int targetTextAlignmentEdit, String strand,
        List<Alignment> alignmentList, List<Integer> endPosList, List<Integer> editNumList) {
        if (strand == "+") {
            targetTextAlignment.setAlignedTarget(targetTextAlignment.getAlignedTargetBuilder().reverse().toString());
            targetTextAlignment.setAlignedTargetBulider(null);
            targetTextAlignment.setAlignedText(targetTextAlignment.getAlignedTextBuilder().reverse().toString());
            targetTextAlignment.setAlignedTextBulider(null);
        }
        else {
            targetTextAlignment.setAlignedTarget(complement(targetTextAlignment.getAlignedTargetBuilder().toString()));
            targetTextAlignment.setAlignedTargetBulider(null);
            targetTextAlignment.setAlignedText(complement(targetTextAlignment.getAlignedTextBuilder().toString()));
            targetTextAlignment.setAlignedTextBulider(null);
        }
        alignmentList.add(targetTextAlignment);
        endPosList.add(targetTextAlignmentPos);
        editNumList.add(targetTextAlignmentEdit);
    }
    
    private static Alignment smithWatermanPostprocessing(
        int[][] M, int m, String strand, String target, String pam, String text, Boolean allowNsInText, Boolean chooseBestInWindow,
        List<Alignment> alignmentList, List<Integer> endPosList, List<Integer> editNumList,
        int col, StringBuilder[] targetTextEmptyAlign, Alignment targetTextAlignment, int targetTextAlignmentPos,
        int targetTextAlignmentEdit) {
        Alignment targetTextAlignmentTemp = find_alignment(M, pam, strand,
            allowNsInText, M[m][col], target, text, m, col, targetTextEmptyAlign, 0, 0);
        if (targetTextAlignmentTemp != null) {
            if (targetTextAlignment == null) {
                return targetTextAlignmentTemp;
            }
            if (((col - targetTextAlignmentPos) > SITE_WINDOW_SIZE)  || !chooseBestInWindow) {
                smithWatermanProcessAlignment(
                    targetTextAlignment, targetTextAlignmentPos, targetTextAlignmentEdit, strand, alignmentList, endPosList, editNumList);

                // return the next possible alignment
                return targetTextAlignmentTemp;
            }
            int targetTextAlignmentTempMisBulge = targetTextAlignmentTemp.getMismatches() + targetTextAlignmentTemp.getBulges();
            int targetTextAlignmentMisBulge = targetTextAlignment.getMismatches() + targetTextAlignment.getBulges();

            if (targetTextAlignmentTempMisBulge < targetTextAlignmentMisBulge) {
                    return targetTextAlignmentTemp;
            }
            if ((targetTextAlignmentTempMisBulge == targetTextAlignmentMisBulge) &&
                (targetTextAlignmentTemp.getBulges() < targetTextAlignment.getBulges())) {
                    return targetTextAlignmentTemp;
                }
        }
        return null;
    }

    private static void smithWatermanLastRowFill(
        int[][] M, int m, int n,  String strand, String target, String pam, String text, int maxEdits,
        int offsetSize, Boolean allowNsInText, Boolean postprocessing, Boolean chooseBestInWindow,
        List<Alignment> alignmentList, List<Integer> endPosList, List<Integer> editNumList) {
        
        Alignment targetTextAlignment = null;
        int targetTextAlignmentPos = 0;
        int targetTextAlignmentEdit = 0;

        // create a empty String Builder for the alignments
        StringBuilder[] targetTextEmptyAlign = null;
        if (postprocessing) {
            targetTextEmptyAlign = new StringBuilder[2];
            targetTextEmptyAlign[0] = new StringBuilder();
            targetTextEmptyAlign[1] = new StringBuilder();
        }

        if (postprocessing) {
            for (int col = 1; col <= n; col++) {
                smithWatermanLastCellLogic(M, target,text, m, col, allowNsInText);
                if ((M[m][col] <= maxEdits) && (col >= offsetSize + 1)) {
                    Alignment targetTextAlignmentTemp = smithWatermanPostprocessing(
                        M, m, strand, target, pam, text, allowNsInText, chooseBestInWindow, alignmentList, endPosList, editNumList,
                        col, targetTextEmptyAlign, targetTextAlignment, targetTextAlignmentPos, targetTextAlignmentEdit);
                    if (targetTextAlignmentTemp != null) {
                        // update the alignment with the next possible alignment
                        targetTextAlignment = targetTextAlignmentTemp;
                        targetTextAlignmentPos = col;
                        targetTextAlignmentEdit = M[m][col];
                    }
                }
            }
            if (targetTextAlignment != null) {
                // TODO: If I'm to solve the window blocks problem (meaning not assuming that each block window is stand alone)
                // then I should provide different soltion for this
                smithWatermanProcessAlignment(
                    targetTextAlignment, targetTextAlignmentPos, targetTextAlignmentEdit, strand, alignmentList, endPosList, editNumList);
            }
        }
        else {
            for (int col = 1; col <= n; col++) {
                smithWatermanLastCellLogic(M, target,text, m, col, allowNsInText);
                if ((M[m][col] <= maxEdits) && (col >= offsetSize + 1)) {
                    endPosList.add(col);
                    editNumList.add(M[m][col]);
                }
            }
        }
    }
    
    public static OffTargetData smithWaterman(
        String strand, String target, String pam, String text, int maxEdits,
        int offsetSize, Boolean allowNsInText, Boolean postprocessing, Boolean chooseBestInWindow) {
        int m = target.length();
        int n = text.length();
        int[][] M = new int[m + 1][n + 1];
        
        // fill the first col
        for (int row = 1; row <= m; row++) {
            M[row][0] = row;
        }
        smithWatermanRowsFill(M, m, n,  target, text, allowNsInText);

        // while filling the last row, we will create the off-target list with the unique end positions
        List<Integer> endPosList = new ArrayList<Integer>();
        List<Integer> editNumList = new ArrayList<Integer>();
        List<Alignment> alignmentList = null;
        if (postprocessing) {
            alignmentList = new ArrayList<Alignment>();
            
        }
        
        smithWatermanLastRowFill(M, m, n, strand, target, pam, text, maxEdits, offsetSize,
           allowNsInText, postprocessing, chooseBestInWindow, alignmentList, endPosList, editNumList);

        
        // convert the lists into arrays
        int[] endPosArr = endPosList.stream().mapToInt(i->i).toArray();
        int[] editNumArr = editNumList.stream().mapToInt(i->i).toArray();
        Alignment[] alignmentArr = null;
        if (postprocessing) {
            alignmentArr = new Alignment[alignmentList.size()];
            alignmentArr = alignmentList.toArray(alignmentArr);
        }

        OffTargetData offTargetData = new OffTargetData(endPosArr.length, endPosArr, editNumArr, alignmentArr);

        return offTargetData;
    }

    private static File[] getChrFiles(String fastaFilePath) throws IOException{
        Path filePath = Path.of(fastaFilePath);
        String fileName = filePath.getFileName().toString();
        String splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
            "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split/";
        File dir = new File(splitFolderPathStr);
        File[] files = dir.listFiles();

        return files;
    }

    public static String reverseComplement(String seq) {
        StringBuilder rcSeq = new StringBuilder();
    
        for (int i = seq.length() - 1; i >= 0; i--) {
            char c = seq.charAt(i);
    
            switch (c) {
                case 'A':
                    rcSeq.append('T');
                    break;
                case 'C':
                    rcSeq.append('G');
                    break;
                case 'G':
                    rcSeq.append('C');
                    break;
                case 'T':
                    rcSeq.append('A');
                    break;
                case 'N':
                    rcSeq.append('N');
                    break;
                case '-':
                    rcSeq.append('-');
                    break;
                default:
                    break;
            }
        }
    
        return rcSeq.toString();
    }

    public static String complement(String seq) {
        StringBuilder cSeq = new StringBuilder();
    
        for (int i = 0; i < seq.length(); i++) {
            char c = seq.charAt(i);
    
            switch (c) {
                case 'A':
                    cSeq.append('T');
                    break;
                case 'C':
                    cSeq.append('G');
                    break;
                case 'G':
                    cSeq.append('C');
                    break;
                case 'T':
                    cSeq.append('A');
                    break;
                case 'N':
                    cSeq.append('N');
                    break;
                case '-':
                    cSeq.append('-');
                    break;
                default:
                    break;
            }
        }

        return cSeq.toString();
    }


    public static class SmithWatermanProcessFileHandler implements Runnable {
        private FileWriter writer;
        private String text;
        private String lastText;
        private int targetLen;
        private Map<String, String> strandToTarget;
        private Map<String, String> strandToPam;
        private Boolean allowNsInText;
        private Boolean postprocessing;
        private Boolean chooseBestInWindow;
        private long positionShift;
        private String chr;
        
        
        public SmithWatermanProcessFileHandler(
            FileWriter writer, String text, String lastText, int targetLen, Map<String, String> strandToTarget, Map<String, String> strandToPam,
            Boolean allowNsInText, Boolean postprocessing, Boolean chooseBestInWindow, long positionShift, String chr) {
            this.writer = writer;
            this.text = text;
            this.lastText = lastText;
            this.targetLen = targetLen;
            this.strandToTarget = strandToTarget;
            this.strandToPam = strandToPam;
            this.allowNsInText = allowNsInText;
            this.postprocessing = postprocessing;
            this.chooseBestInWindow = chooseBestInWindow;
            this.positionShift = positionShift;
            this.chr = chr;
        }
        
        @Override
        public void run() {
            int offsetSize = 0;
            if (lastText != null){
                offsetSize = targetLen + MAX_EDITS + 1; // targetLen + MAX_EDITS should be fine as well
                text = lastText.substring(lastText.length() - offsetSize) + text;
            }

            for (String strand : strandToTarget.keySet()) {
                OffTargetData offTargetData = smithWaterman(
                    strand, strandToTarget.get(strand), strandToPam.get(strand), text, MAX_EDITS, offsetSize, allowNsInText, postprocessing, chooseBestInWindow);
                // write to csv
                int size = offTargetData.getSize();
                synchronized(writer){
                try {
                    for (int i = 0; i < size; i++) {
                        // TODO: find solution for the text.replace("\n", "").replace("\r", "")
                        // write chromosome and starnd
                        writer.append(chr + "," + strand + ",");
                        // write end position and the site with flank
                        int relativePos = offTargetData.getEndPosArr()[i];
                        writer.append(String.valueOf(relativePos + positionShift - offsetSize));
                        writer.append(",");
                        writer.append(text.substring(Math.max(relativePos - targetLen - MAX_EDITS, 0), relativePos).replace("\n", "").replace("\r", ""));
                        writer.append(",");
                        // write edit size
                        writer.append(String.valueOf(offTargetData.getEditNumArr()[i]));
                        if (postprocessing) {
                            writer.append(",");
                            // write aligned target
                            Alignment alignment = offTargetData.getAlignmentArr()[i];
                            writer.append(alignment.getAlignedTarget());
                            writer.append(",");
                            // write aligned text
                            writer.append(alignment.getAlignedText().replace("\n", "").replace("\r", ""));
                            writer.append(",");
                            // write aligment mismatches
                            writer.append(String.valueOf(alignment.getMismatches()));
                            writer.append(",");
                            // write aligment bulges
                            writer.append(String.valueOf(alignment.getBulges())); 
                        }
                        writer.append("\n");
                    }
                    writer.flush(); 
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }
            }
        }
    }
    

    public static void smithWatermanProcessFile(
            File file, Map<String, String> strandToTarget, Map<String, String> strandToPam,
            FileWriter writer, int bufferSize, int targetLen,
            Boolean allowNsInText, Boolean postprocessing, Boolean chooseBestInWindow) throws IOException{
        Path filePath = Path.of(file.getAbsolutePath());
        // Get the chr name
        String chr = file.getName();
        chr = chr.substring(0, chr.lastIndexOf('.'));
        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);
        List<Future<?>> futures = new ArrayList<>();
        InputStream in = Files.newInputStream(filePath);
        long positionShift = 0;
        byte[] buffer = new byte[bufferSize];
        int numRead;
        // int iter = 1;
        String lastText = null;
        while ((numRead = in.readNBytes(buffer, 0, bufferSize)) != 0) {
            if (numRead == -1) {
                break;
            }
            String text = new String(buffer, 0, numRead);
            SmithWatermanProcessFileHandler handler = new SmithWatermanProcessFileHandler(
                writer, text, lastText, targetLen, strandToTarget, strandToPam, allowNsInText,
                postprocessing, chooseBestInWindow, positionShift, chr);
            Future<?> future = executorService.submit(handler);
            futures.add(future);
            
            // System.out.println(iter);
            // iter ++;
            positionShift += bufferSize;
            lastText = text;
            }
            try {
                for (Future<?> future : futures) {
                    future.get(); // wait for all tasks to complete
                }
            }
            catch (Exception e){
                e.printStackTrace();
            }
            // Close input stream
            in.close();
            executorService.shutdown();
    }

    public static void main(String[] args) {
        // String[] args1 = {"genomes/hg38_only_chrs.fa", "sgRNAs.txt", "output_prefix", MaxEdits, MaxMismatchesWithoutBulges, MaxMismatchesWithBulges, MaxBulges, NUM_THREADS, chooseBestInWindow, best_in_window_size, PAM, ALLOW_PAM_EDITS};
        // args = args1;

        // Start measuring execution time
        // String PAM = "NGG";
        Instant start = Instant.now();
        
        // Some constants that should be provided as args
        Boolean allowNsInText = false;
        Boolean postprocessing = true;
        String fastaFilePath = args[0];

        String filePath = args[1];
        Boolean isFileInput = true;
        
        // Set the values of the static variables
        SmithWatermanOffTargetSearchAlign.setMaxEdits(Integer.parseInt(args[3]));
        SmithWatermanOffTargetSearchAlign.setMaxMismatchesWithoutBulges(Integer.parseInt(args[4]));
        SmithWatermanOffTargetSearchAlign.setMaxMismatchesWithBulges(Integer.parseInt(args[5]));
        SmithWatermanOffTargetSearchAlign.setMaxBulges(Integer.parseInt(args[6]));
        SmithWatermanOffTargetSearchAlign.setNumThreads(Integer.parseInt(args[7]));
        Boolean chooseBestInWindow = args[8].equals("true");
        SmithWatermanOffTargetSearchAlign.setSiteWindowSize(Integer.parseInt(args[9]));
        String PAM = args[10]; // "NGG";
        SmithWatermanOffTargetSearchAlign.setAllowPamEdits(args[11].equals("true"));


        Map<String, String> strandToPam = new HashMap<>();
        strandToPam.put("+", PAM);
        strandToPam.put("-", reverseComplement(PAM));

        try {
            if (filePath.indexOf('.') == -1){
                FileWriter writerTemp = new FileWriter("target_temp.txt");
                writerTemp.append(filePath + "\n");
                writerTemp.close();
                filePath = "target_temp.txt";
                isFileInput = false;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String target;
            while ((target = reader.readLine()) != null) {
                String rcTarget = reverseComplement(target);
                // This is import factor as it define the memory the program consumes
                // In addition, we are limited by the size of the array which function of this size
                int bufferSize = 10000;

                Map<String, String> strandToTarget = new HashMap<>();
                strandToTarget.put("+", target);
                strandToTarget.put("-", rcTarget);
                int targetLen = target.length();
                // Split genome to chrs files (if not already done)
                FastaReader.splitFastaToFiles(fastaFilePath);
                File[] files = getChrFiles(fastaFilePath);

                // Define output file and write col names
                FileWriter writer;
                if (isFileInput){
                    writer = new FileWriter(args[2] + "_" + target + ".csv");
                }
                else{
                    writer = new FileWriter(args[2] + ".csv");
                }

                String[] cols;
                if (postprocessing) {
                    cols = new String[] {CHROMOSOME_STR, STRAND_STR, END_POS_STR, SEQ_STR, EDIT_NUM_STR,
                        ALINGED_TARGET_STR, ALINGED_TEXT_STR, MISMATCHES_NUM_STR, BUGLES_NUM_STR};
                }
                else{
                    cols = new String[] {CHROMOSOME_STR, STRAND_STR, END_POS_STR, SEQ_STR, EDIT_NUM_STR,
                        ALINGED_TARGET_STR, ALINGED_TEXT_STR, MISMATCHES_NUM_STR, BUGLES_NUM_STR};
                }
                for (int colIndex = 0; colIndex < cols.length; colIndex++) {
                    writer.append(cols[colIndex]);
                    if (colIndex < (cols.length - 1)) {
                        writer.append(",");
                    }
                }
                writer.append("\n");

                for (int file_i=0 ; file_i < files.length; file_i++) {
                    smithWatermanProcessFile(
                        files[file_i], strandToTarget, strandToPam, writer, bufferSize, targetLen,
                        allowNsInText, postprocessing, chooseBestInWindow);
                }
                // Close writer
                writer.close();
                // hint ;the garbage collector
                System.gc();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        (new File("target_temp.txt")).delete();
        // report execution time
        Instant end = Instant.now();
        Duration timeElapsed = Duration.between(start, end);
        System.out.println("total Execution time: " + timeElapsed.toMillis() + " milliseconds");
    }
}
