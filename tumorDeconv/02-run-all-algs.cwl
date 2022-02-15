#!/usr/bin/env cwltool
class: Workflow
label: 02-run-all-algs
id: run-all-algos
cwlVersion: v1.2

requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement

inputs:
   signatures:
     type: string[]
   prot-algorithms:
     type: string[]
   inFile:
     type: File

outputs:
   deconvoluted:
     type: File[]
     outputSource:
      - run-all-algs-by-sig/deconvoluted


steps:
   get-all-mat:
      run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/main/signature_matrices/get-signature-matrix.cwl
      #./../proteomicsTumorDeconv/signature_matrices/get-signature-matrix.cwl
      scatter: [sigMatrixName]
      scatterMethod: flat_crossproduct
      in:
        sigMatrixName: signatures
      out:
        [sigMatrix]
   run-all-algs-by-sig:
      run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/main/metrics/run-deconv.cwl
      #./../proteomicsTumorDeconv/metrics/run-deconv.cwl
      scatter: [signature,alg]
      scatterMethod: flat_crossproduct    
      in:
        signature: get-all-mat/sigMatrix
        alg: prot-algorithms
        matrix: inFile
      out:
        [deconvoluted]
