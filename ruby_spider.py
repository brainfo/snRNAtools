#coding=utf-8
from urllib.request import urlopen
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup
import re,os,sys
import time

 from scrapy.spiders import CrawlSpider, Rule, BaseSpider, Spider
 from scrapy.linkextractors.lxmlhtml import LxmlLinkExtractor 
 from scrapy.selector import Selector
 from scrapy.http import HtmlResponse

## this is a test script for crawling geo data ##
#########     2019-08-08 Ruby Jiang     #########

class MySpider(scrapy.Spider):
    name = 'myspider'
    start_urls = ['http://catalog.umassd.edu/content.php?catoid=45&navoid=3554']

    def parse(self, response):
        # get "onclick" java function of every "show more" link
        # and extract parameters supplied to this function with regular expressions
        links = response.xpath("//a/@onclick[contains(.,'showHide')]")
        for link in links:
            args = link.re("'(.+?)'")
            # make our url by putting arguments from page source 
            # into a template of an url
            url = 'http://catalog.umassd.edu/ajax/preview_filter_show_hide_data.php?show_hide={}&cat_oid={}&nav_oid={}&ent_oid={}&type={}&link_text={}'.format(*args)
            yield scrapy.Request(url, self.parse_more) 

    def parse_more(self, response):
        # here you'll get page source with all of the links