#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GERMLINE } from './workflows/germline'

workflow {
    GERMLINE()
}
