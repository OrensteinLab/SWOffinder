package SmithWatermanOffTarget;

public class OffTargetData {
    private int size;
    private int[] endPosArr;
    private int[] editNumArr;
    private Alignment[] alignmentArr;

    public OffTargetData(int size, int[] endPosArr, int[] editNumArr, Alignment[] alignmentArr) {
        this.size = size;
        this.endPosArr = endPosArr;
        this.editNumArr = editNumArr;
        this.alignmentArr = alignmentArr;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }
    
    public int[] getEndPosArr() {
        return endPosArr;
    }

    public void setEndPosArr(int[] endPosArr) {
        this.endPosArr = endPosArr;
    }

    public int[] getEditNumArr() {
        return editNumArr;
    }

    public void setEditNumArr(int[] editNumArr) {
        this.editNumArr = editNumArr;
    }

    public Alignment[] getAlignmentArr() {
        return alignmentArr;
    }

    public void setAlignmentArr(Alignment[] alignmentArr) {
        this.alignmentArr = alignmentArr;
    }
}