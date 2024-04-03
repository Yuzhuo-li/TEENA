#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1 
# Update 12-09: We have set an option to evaluate the effectiveness of region set similarity metrics. For example, we can evaluate differences in behavior between the query file (or a TE repbase file) that removes the gap or promoter regions of the genome sequence between their original file.

import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import euclidean_distances, cosine_similarity
from scipy.spatial.distance import jaccard

# 
def calculate_metrics(file_name, file1, file2):
    with open(file1, 'r') as file:
        doc1 = file.read()
    with open(file2, 'r') as file:
        doc2 = file.read()

    # Convert text to word frequency vector
    vectorizer = CountVectorizer()
    X = vectorizer.fit_transform([doc1, doc2]).toarray()

    # Euclidean distance
    euclidean_dist = euclidean_distances(X[0].reshape(1, -1), X[1].reshape(1, -1))

    # Cosine similarity
    cosine_sim = cosine_similarity(X[0].reshape(1, -1), X[1].reshape(1, -1))

    # Convert to a set to calculate coverage score and jaccard distance
    set1 = set(doc1.split())
    set2 = set(doc2.split())

    # Coverage score
    overlap_score = len(set1 & set2) / min(len(set1), len(set2))

    # Jaccard distance
    jaccard_dist = jaccard(X[0], X[1])

   
    f = open(file_name, 'w')
    f.write(f"Euclidean Distance: {euclidean_dist[0][0]} \n")
    f.write(f"Cosine Similarity: {cosine_sim[0][0]} \n")
    f.write(f"Overlap Score: {overlap_score} \n")
    f.write(f"Jaccard Distance: {jaccard_dist}")
    f.flush()
    f.close()
   
    
