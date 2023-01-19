import math
import selenium
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By 
# options = webdriver.ChromeOptions()
# options.add_experimental_option('excludeSwitches', ['enable-logging'])
# wd = webdriver.Chrome(options=options)
# wd = webdriver.Firefox()
#lN*E$7b&9ZsrRJ4TV96t8
wd = webdriver.Chrome('./chromedriver.exe')


search_numbers = """3.2.1.4
3.2.1.91
3.2.1.21
3.2.1.8
3.2.1.37
3.2.1.78
3.2.1.25
3.2.1.55
3.2.1.99
3.2.1.131"""

final_dict = {"Protein name":[], "EC#":[], "Uniprot":[], 'name':[]}

for search_number in search_numbers.split():
    print(search_number)
    initial_url = f'http://www.cazy.org/search?page=recherche&lang=en&recherche={search_number}&tag=9'
    wd.get(initial_url)
    # elem = wd.find_element("class_name","paragraphe")
    # print(elem.text)
    

    aTagsInLi = wd.find_elements(By.CSS_SELECTOR,'li a')
    element_list = []
    hits = wd.find_element(By.CSS_SELECTOR,'h2')
    hits = math.ceil(int(hits.text.split('hits')[0].split(' ')[-2][1:])/10)
    # print(hits)
    for a in aTagsInLi:
        if str(a.get_attribute('href')).endswith('html'):
            element_list.append(a.get_attribute('href'))
    for i in range(hits-1):

        initial_2_url = f'{initial_url}&debut_knoact={(i+1)*10}#pagination_knoact'
        wd.get(initial_2_url)
        aTagsInLi = wd.find_elements(By.CSS_SELECTOR,'li a')
        for a in aTagsInLi:
            if str(a.get_attribute('href')).endswith('html'):
                element_list.append(a.get_attribute('href'))

    print(element_list)

    big_dict = {}
    for elem in element_list:
        name = elem.split('/')[-1].split('.html')[0]
        elem = f"{elem.split('.html')[0]}_characterized.html"
        print(elem)
        big_dict[name] = final_dict
        wd.get(elem)
        tab = wd.find_elements(By.TAG_NAME,'table')
        # print(tab)
        for i in tab[1:]:
            rows = i.find_elements(By.TAG_NAME,'tr')
            for row in rows:
                tds = row.find_elements(By.TAG_NAME,'td')
                if len(tds)>4 and tds[1].text == search_number:
                    final_dict['Protein name'].append(tds[0].text)
                    final_dict['EC#'].append(tds[1].text)
                    final_dict['Uniprot'].append(tds[5].text)
                    final_dict['name'].append(name)

df = pd.DataFrame.from_dict(final_dict)
# print(big_dict)
df
df.to_csv('fin.csv')