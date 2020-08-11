import os
import shutil
import nltk
from nltk.corpus import stopwords
from nltk import pos_tag
import random
from nltk.stem.lancaster import LancasterStemmer
from nltk import pos_tag 
import sys
import numpy as np
from gensim.models import Word2Vec
question_per_exam = [7, 7, 7, 7, 4, 7, 7, 7, 7, 7, 10, 10]
max_user = 24
answers_path = "answers"
scores_path = "scores"
keywords_path = "keywords"

keep_tagger = set([
'JJ','JJR','JJS', 
'CD',
'NN',
'NNS','NNP','NNPS',
'PDT',
'SYM',
'RB','RBR','RBS'
'VB','VBD',
'VBG','VBN','VBP','VBZ'
	])

fs = os.listdir(answers_path)


STOP_WORDS = stopwords.words('english')
STOP_WORDS = set(STOP_WORDS)
STOP_WORDS.remove('both')
STOP_WORDS.remove('all')
STOP_WORDS.add('br')
# make dir
# question_dir_name = ['q'+str(i+1) for i in range(12)]
# for qd in question_dir_name:
# 	adirpath = os.path.join(answers_path,qd)
# 	os.mkdir(adirpath)
# 	sdirpath = os.path.join(scores_path,qd)
# 	os.mkdir(sdirpath)
# 	kdirpath = os.path.join(keywords_path,qd)
# 	os.mkdir(kdirpath)


# move file
# for base_dir in [answers_path, scores_path]:
# 	fs = os.listdir(base_dir)
# 	fs = [f for f in fs if not f.startswith('q')]
# 	question_per_exam_counter = [7]
# 	for i in range(1,len(question_per_exam)):
# 		question_per_exam_counter.append(question_per_exam_counter[i-1]+question_per_exam[i])
# 	# print(question_per_exam_counter)

# 	for f in fs:
# 		fp = os.path.join(base_dir,f)
# 		fid = int(f)+1
# 		qid = -1
# 		for i in range(len(question_per_exam_counter)):
# 			if fid <= question_per_exam_counter[i]:
# 				qid = i + 1
# 				break
# 		print(fid, qid)
# 		qp = "q"+str(qid)
# 		qp = os.path.join(base_dir,qp,str(fid-1))
# 		shutil.move(fp,qp)
	# break


# read answers

pre_trained_vec_path = "glove.6B.50d.txt"

tokenizer = nltk.RegexpTokenizer(r"\w+").tokenize
st = LancasterStemmer()
# question_dir_name = ['q'+str(i+1) for i in range(12)]
answer_dir_path =  answers_path
keywords_dir_path =  keywords_path

primitive_pre_trained_dic = {}

with open(pre_trained_vec_path) as f:
	line = f.readline()
	while line:
		line = line.split()
		w = line[0]
		vec = ' '.join(line[1:])
		primitive_pre_trained_dic[w]=vec 
		line = f.readline()


# filter dictionary
# if not os.path.exists(model_path):
files = os.listdir(answer_dir_path)
keywords_set = set()
for fname in files:
	fpath = os.path.join(answer_dir_path,fname)
	with open(fpath,'rb') as f:
		lines = f.readlines()
		lines = lines[:max_user]
		for line in lines:
			line = line.decode('latin-1')
			line = line.strip()
			keywords = tokenizer(line)

			word_and_tagger = pos_tag(keywords)
			# print(word_and_tagger)
			keywords = [e[0] for e in word_and_tagger if e[1] in keep_tagger]
			keywords = [w.lower() for w in keywords]
			keywords = [w for w in keywords if w not in STOP_WORDS]
			# keywords = [st.stem(w) for w in keywords]
			keywords = [w for w in keywords if w in primitive_pre_trained_dic]
			keywords_set.update(set(keywords))


pre_trained_dic = {}
for k,v in primitive_pre_trained_dic.items():
	if k in keywords_set:
		pre_trained_dic[k] = v 

allwords = list(keywords_set)


# train word2vec
# get corpus
# model_path = "word2vec_model"
# # if not os.path.exists(model_path):
# corpus = []
# files = os.listdir(answer_dir_path)
# for fname in files:
# 	fpath = os.path.join(answer_dir_path,fname)
# 	with open(fpath,'rb') as f:
# 		lines = f.readlines()
# 		lines = lines[:max_user]
# 		for line in lines:
# 			line = line.decode('latin-1')
# 			line = line.strip()
# 			keywords = tokenizer(line)

# 			word_and_tagger = pos_tag(keywords)
# 			# print(word_and_tagger)
# 			keywords = [e[0] for e in word_and_tagger if e[1] in keep_tagger]
# 			keywords = [w.lower() for w in keywords]
# 			keywords = [w for w in keywords if w not in STOP_WORDS]
# 			keywords = [st.stem(w) for w in keywords]
# 			if len(keywords)!=0:
# 				corpus.append(keywords)
				
# model = Word2Vec(corpus, size = 100, window = 2, min_count = 1, workers = 4)
# model.save(model_path)
# model = Word2Vec.load(model_path)	


# get statistics
c = 0
key_num = []
ans_len = []
#get maximum length
files = os.listdir(answer_dir_path)
for fname in files:
	fpath = os.path.join(answer_dir_path,fname)
	with open(fpath,'rb') as f:
		lines = f.readlines()
		lines = lines[:max_user]
		for line in lines:
			line = line.decode('latin-1')
			ans_len.append(len(line.split()))

print("maximum answer length", max(ans_len))


#extract keywords

files = os.listdir(answer_dir_path)
for fname in files:
	fpath = os.path.join(answer_dir_path,fname)
	keyword_content = []
	with open(fpath,'rb') as f:
		lines = f.readlines()
		lines = lines[:max_user]
		for line in lines:
			line = line.decode('latin-1')
			line = line.strip()
			keywords = tokenizer(line)

			word_and_tagger = pos_tag(keywords)
			keywords = [e[0] for e in word_and_tagger if e[1] in keep_tagger]
			keywords = [w.lower() for w in keywords]
			keywords = [w for w in keywords if w not in STOP_WORDS]
			# keywords = [st.stem(w) for w in keywords]
			keywords = [w for w in keywords if w in primitive_pre_trained_dic]

			if len(keywords) == 0:
				print(line)
				c+=1
				# keywords = ['None']
				keywords = random.sample(allwords,1)
			keywords = list(set(keywords))
			key_num.append(len(keywords))
			keyword_content.append(' '.join(keywords) + '\n')
	kpath = os.path.join(keywords_dir_path, fname)
	with open(kpath,'w') as f:
		f.writelines(keyword_content)
print("keywords statistics, min, max, avg:",min(key_num),max(key_num),sum(key_num)/(len(key_num)))
print("neglect answer number:",c)

print("length of answer_dictionary",len(keywords_set))


# model.wv.save("pretrained_vec")

word_vec_path = "pretrained_vec"
lines = []
with open(word_vec_path,'w') as f:
	for k,v in pre_trained_dic.items():
		k = k.lower()
		# if k in STOP_WORDS:
		# 	continue
		v = v.split()
		v = np.array([eval(e) for e in v])
		v = v/np.linalg.norm(v)
		v = [str(e) for e in v]
		lines.append(k + ' ' + ' '.join(v) + '\n')
	f.writelines(lines)
print("finish write dictioanry")
			
