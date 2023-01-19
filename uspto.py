import math
import selenium
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import time 
from bs4 import BeautifulSoup as bs
# options = webdriver.ChromeOptions()
# options.add_experimental_option('excludeSwitches', ['enable-logging'])
# wd = webdriver.Chrome(options=options)
# wd = webdriver.Firefox()
#lN*E$7b&9ZsrRJ4TV96t8
wd = webdriver.Chrome('./chromedriver.exe')
doc_number = '10596246'
url = f"https://patft.uspto.gov/netacgi/nph-Parser?TERM1={doc_number}&Sect1=PTO1&Sect2=HITOFF&d=PALL&p=1&u=%2Fnetahtml%2FPTO%2Fsrchnum.htm&r=0&f=S&l=50"
wd.get(url)
# tab = wd.find_element(By.XPATH,"//td[contains(text(),'Current CPC Class:')]")
# tab = wd.find_element(By.TAG_NAME,'p')
# print(tab.text)
# # print(tab)
# for i in tab[1:]:
#     rows = i.find_elements(By.TAG_NAME,'tr')
#     for row in rows:
#         print(row.text)
time.sleep(3)
# tab = wd.find_element(By.XPATH,"//td[contains(text(),'Current CPC Class:')]")
l = wd.find_element(By.CSS_SELECTOR,"html")
text = l.get_attribute('innerHTML')
text = str(text)
ind = text.find("Current CPC Class")

# print(text[ind+20:ind+1000])
text = text[text.find('Current CPC Class:'):text.find('Current CPC Class:')+1000]
soup = bs(text, 'html.parser')
ds = soup.findall('td')
print(ds)
# print("HTML code of element: " + l.get_attribute('innerHTML'))