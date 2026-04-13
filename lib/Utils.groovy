class Utils {

    /**
     * Look up a genome attribute from params.genomes (loaded via igenomes.config).
     */
    static String getGenomeAttr(def params, String attr) {
        if (!params.genome) return null
        if (!params.genomes?.containsKey(params.genome)) {
            throw new Exception("Unknown --genome '${params.genome}'. Add it to conf/igenomes.config")
        }
        def g = params.genomes[params.genome]
        if (!g.containsKey(attr)) {
            throw new Exception("Genome '${params.genome}' does not define '${attr}' in conf/igenomes.config")
        }
        return g[attr]
    }

    /**
     * Parse sample → patient mappings from a CSV samplesheet.
     * Returns [sampleIds: List<String>, sample2patient: Map<String,String>].
     */
    static Map readSampleSheetInfo(def csvPathObj) {
        def csvPath = csvPathObj.toString()
        def f = new File(csvPath)
        if (!f.exists()) throw new Exception("Samplesheet not found: ${csvPath}")

        def lines = f.readLines().findAll { it?.trim() }
        if (lines.size() < 2) throw new Exception("Samplesheet has no data rows: ${csvPath}")

        def header     = lines[0].split(',', -1)*.trim()
        def idxSample  = header.indexOf('sample')
        def idxPatient = header.indexOf('patient')

        if (idxSample < 0) throw new Exception("Samplesheet must have a 'sample' column: ${csvPath}")

        def sample2patient = [:]
        lines.tail().each { ln ->
            def cols = ln.split(',', -1)
            if (cols.size() > idxSample) {
                def sid = cols[idxSample]?.trim()
                if (sid && !sample2patient.containsKey(sid)) {
                    def pat = (idxPatient >= 0 && idxPatient < cols.size()) ? cols[idxPatient]?.trim() : ''
                    sample2patient[sid] = pat ?: sid
                }
            }
        }

        if (sample2patient.isEmpty()) throw new Exception("No sample IDs parsed from: ${csvPath}")

        return [
            sampleIds      : sample2patient.keySet().toList().sort(),
            sample2patient : sample2patient
        ]
    }

    /**
     * Returns "bam" if the samplesheet has a 'bam' column, "fastq" otherwise.
     */
    static String detectSampleSheetType(def csvPathObj) {
        def csvPath = csvPathObj.toString()
        def f = new File(csvPath)
        if (!f.exists()) throw new Exception("Samplesheet not found: ${csvPath}")
        def lines = f.readLines().findAll { it?.trim() }
        if (lines.isEmpty()) throw new Exception("Samplesheet is empty: ${csvPath}")
        def header = lines[0].split(',', -1)*.trim()*.toLowerCase()
        return header.contains('bam') ? 'bam' : 'fastq'
    }

    /** Returns sample IDs missing a .bam or .bam.bai in bamDir. */
    static List findMissingBamBai(List sampleIds, String bamDir) {
        sampleIds.findAll { sid ->
            !new File("${bamDir}/${sid}.bam").exists() ||
            !new File("${bamDir}/${sid}.bam.bai").exists()
        }
    }

    /** Returns sample IDs missing a .metrics.txt in bamDir. */
    static List findMissingMetrics(List sampleIds, String bamDir) {
        sampleIds.findAll { sid -> !new File("${bamDir}/${sid}.metrics.txt").exists() }
    }
}
