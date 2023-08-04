package SmithWatermanOffTarget;

public class Alignment {
    private int mismatches;
    private int bulges;
    private StringBuilder alignedTargetBulider;
    private StringBuilder alignedTextBulider;
    private String alignedTarget = null;
    private String alignedText = null;

    public Alignment(int mismatches, int bulges, StringBuilder alignedTargetBulider, StringBuilder alignedTextBulider) {
        this.mismatches = mismatches;
        this.bulges = bulges;
        this.alignedTargetBulider = alignedTargetBulider;
        this.alignedTextBulider = alignedTextBulider;
    }

    public int getMismatches() {
        return mismatches;
    }

    public void setMismatches(int mismatches) {
        this.mismatches = mismatches;
    }

    public int getBulges() {
        return bulges;
    }

    public void setBulges(int bulges) {
        this.bulges = bulges;
    }
    
    public StringBuilder getAlignedTargetBuilder() {
        return alignedTargetBulider;
    }

    public void setAlignedTargetBulider(StringBuilder alignedTargetBulider) {
        this.alignedTargetBulider = alignedTargetBulider;
    }

    public StringBuilder getAlignedTextBuilder() {
        return alignedTextBulider;
    }

    public void setAlignedTextBulider(StringBuilder alignedTextBulider) {
        this.alignedTextBulider = alignedTextBulider;
    }

    public String getAlignedTarget() {
        return alignedTarget;
    }

    public void setAlignedTarget(String alignedTarget) {
        this.alignedTarget = alignedTarget;
    }

    public String getAlignedText() {
        return alignedText;
    }

    public void setAlignedText(String alignedText) {
        this.alignedText = alignedText;
    }

}
