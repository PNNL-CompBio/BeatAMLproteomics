#!/usr/bin/env cwltool
class: Workflow
label: 02-run-all-algs
id: run-all-algos
cwlVersion: v1.2

requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement

inputs:
   signatures:
     type: string[]
   prot-algorithms:
     type: string[]
   protFile:
     type: File
   rnaFile:
     type: File

outputs:
   deconvoluted:
     type: File
     outputSource:
      - run-best-algs-by-sig/deconvoluted

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
   get-all-cors:
      run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/main/metrics/correlations/deconv-corr-cwl-tool.cwl
      scatter: [signature,mrnaAlg,protAlg]
      scatterMethod: flat_crossproduct
      in:
        transcriptomics: rnaFile
        proteomics: protFile
        cancerType:
           valueFrom: "AML"
        mrnaAlg: prot-algorithms
        protAlg: prot-algorithms
        signature: get-all-mat/sigMatrix
        sampleType:
           valueFrom: "all"
      out:
        [corr]
   get-best-cor:
       run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/localdata/metrics/correlations/best-deconv-cor-tool.cwl
       in:
         corFiles: get-all-cors/corr
       out:
         value
   get-best-mat:
       run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/main/signature_matrices/get-signature-matrix.cwl
       in:
          sigMatrixName: get-best-cor/value[1]
       out:
          [sigMatrix]
   run-best-algs-by-sig:
      run: https://raw.githubusercontent.com/PNNL-CompBio/proteomicsTumorDeconv/main/metrics/run-deconv.cwl
      in:
        signature: get-best-mat/sigMatrix
        alg: get-best-cor/value[2]
        matrix: protFile
      out:
        [deconvoluted]
