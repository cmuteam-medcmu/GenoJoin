import faiss
import numpy as np
from time import time
# from nltk.tokenize import word_tokenize
from gensim.models import Word2Vec

# Example sentences
# sentences_train = [
#     "CHR1 10000001 C A.",
#     "CHR1 10000001 C G.",
#     "CHR1 10000001 C T.",
#     "CHR1 10000002 T C.",
#     "CHR1 10000002 T T.",
#     "CHR1 10000005 A AG.",
#     "CHR1 10000012 C -.",
# ]

class WordDatabase:
    def __init__(self, sentences_train):
        self.sentences_train = sentences_train

    @classmethod
    def sentence_vector(cls, tokens, model):
        vecs = [model.wv[word] for word in tokens if word in model.wv]
        if len(vecs) == 0:
            return np.zeros(model.vector_size)
        return np.mean(vecs, axis=0)

    def CreateDB(self, threads):
        start_time = time()

        # Tokenize
        # tokenized_sentences = [word_tokenize(s) for s in self.sentences_train]
        tokenized_sentences = [s.split() for s in self.sentences_train]

        # Train Word2Vec
        self.model = Word2Vec(tokenized_sentences, vector_size=300, window=5, min_count=1, workers=threads, seed=42)


        sentence_vectors = np.array([self.sentence_vector(t, self.model) for t in tokenized_sentences], dtype='float32')

        d = self.model.vector_size # vector dimension
        
        faiss.normalize_L2(sentence_vectors)  
        self.index = faiss.IndexFlatIP(d)  
        self.index.add(sentence_vectors)

        end_time = time() - start_time
        print(f'>> Create Database: {end_time:.2f} sec')


    def QueryDB(self, sentences_query=None):
        start_time = time()
        # sentences_query = [
        #     'chr16 994022 T A', 
        #     'chr16 994022 T T',
        #     'chr16 994022 T CA',
        #     'chr16 994022 T T',
        # ]
        # n = len(sentences_query)
        # by = 10
        # cy = (n//by)+1
        # start = 0
        # end = by
        # for i in range(cy):
        start_ti = time()
        tokenized_query = [s.split() for s in sentences_query[0:10]]
        vectors_query = np.array([self.sentence_vector(t, self.model) for t in tokenized_query], dtype='float32')
        faiss.normalize_L2(vectors_query)

            # --- Search top 3 similar sentences in the training set
        distances, indices = self.index.search(vectors_query, 1)

            # start = end 
            # if i+1 != cy:
            #     end = by * (i+2)
            # else:
            #     end = n
            # print(start, end)

            # print(f'\n>> Cyc {i+1}: {(time() - start_ti):.2f} sec')
        # print(indices)
        # print(distances)

        for i, (idxs, sims) in enumerate(zip(indices, distances)):
            # print(f"\nQuery sentence {i}: '{sentences_query[i]}'")
            for rank, (j, sim) in enumerate(zip(idxs, sims)):
                if sentences_query[i] != self.sentences_train[j]:
                    print(f"\nQuery sentence '{sentences_query[i]}' => False")
                else:
                    print(f"\nQuery sentence '{sentences_query[i]}' => True")
                print(f"  sentence {j} ('{self.sentences_train[j]}'), similarity={sim:.2f}\n")

        end_time = time() - start_time
        print(f'>> Query Database: {end_time:.2f} sec')
        # return model, index