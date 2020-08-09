import os
import shutil
import nltk
from nltk.corpus import stopwords
from nltk import pos_tag
import random
from nltk.stem import WordNetLemmatizer 
import sys
from gensim.models import Word2Vec
question_per_exam = [7, 7, 7, 7, 4, 7, 7, 7, 7, 7, 10, 10]
max_user = 24
answers_path = "answers"
scores_path = "scores"
keywords_path = "keywords"

fs = os.listdir(answers_path)
fs = [f for f in fs if not f.startswith('q')]

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

# load dictionary 
pretrained_vector_path = "glove.twitter.27B.50d.txt" 
vec_dictionary = set()
pretrained_vec = dict()
with open(pretrained_vector_path) as f:
	l = f.readline()
	c = 0
	while l:
		# if c>=20:
		# 	sys.exit(-1)
		# print(l)
		l = l.split()
		w = l[0]
		vec = ' '.join(l[1:])
		pretrained_vec[w] = vec
		vec_dictionary.add(w)
		l = f.readline()
		c+=1
# read answers

tokenizer = nltk.tokenize.WordPunctTokenizer().tokenize
lemmatizer = WordNetLemmatizer().lemmatize
question_dir_name = ['q'+str(i+1) for i in range(12)]
answer_dir_path = [os.path.join(answers_path, f) for f in question_dir_name]
keywords_dir_path = [os.path.join(keywords_path, f) for f in question_dir_name]


# get dictioanry
answer_dictionary = set()
for d in answer_dir_path:
	files = os.listdir(d)
	for fname in files:
		fpath = os.path.join(d,fname)
		with open(fpath,'rb') as f:
			lines = f.readlines()
			lines = lines[:max_user]
			for line in lines:
				line = line.decode('latin-1')
				line = line.strip()
				keywords = tokenizer(line)
				# keywords = [lemmatizer(w) for w in keywords]
				keywords = [w for w in keywords if w not in stopwords.words('english')]
				answer_dictionary.update(set(keywords))
# remove words not in vec_dictionary
print("length of answer_dictionary",len(answer_dictionary))
tmp_answer_dictionary = set()
for w in answer_dictionary:
	if w in vec_dictionary:
		tmp_answer_dictionary.add(w)
answer_dictionary = tmp_answer_dictionary
print("after remove words not in vec, length of answer_dictionary",len(answer_dictionary))


# get statistics
c = 0
key_num = []
keywords_list = list(answer_dictionary)
for d, k in zip(answer_dir_path, keywords_dir_path):
	files = os.listdir(d)
	for fname in files:
		fpath = os.path.join(d,fname)
		keyword_content = []
		with open(fpath,'rb') as f:
			lines = f.readlines()
			lines = lines[:max_user]
			for line in lines:
				line = line.decode('latin-1')
				line = line.strip()
				keywords = tokenizer(line)
				# keywords = [lemmatizer(w) for w in keywords]
				keywords = [w for w in keywords if w not in stopwords.words('english')]
				tmp_keywords = keywords.copy()
				keywords = [w for w in keywords if w not in answer_dictionary]
				if len(keywords) == 0:
					c+=1
					print(tmp_keywords)
					keywords = random.sample(keywords_list,1)
				key_num.append(len(keywords))
				keyword_content.append(' '.join(keywords) + '\n')
		kpath = os.path.join(k, fname)
		with open(kpath,'w') as f:
			f.writelines(keyword_content)
print("keywords statistics, min, max, avg:",min(key_num),max(key_num),sum(key_num)/(len(key_num)))
print("neglect answer number:",c)

word_vec_path = "pretrained_vec"
lines = []
with open(word_vec_path,'w') as f:
	for k,v in pretrained_vec.items():
		if k not in answer_dictionary:
			continue
		else:
			lines.append(k + ' ' + v + '\n')
	f.writelines(lines)

			


