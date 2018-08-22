#!/bin/bash

tag=$1

sed -i "" -e "s#datasources/#/datasources/#" common_variant_filter.py

docker build -t vanallenlab/common_variant_filter:${tag} .

sed -i "" -e "s#/datasources/#datasources/#" common_variant_filter.py

docker push vanallenlab/common_variant_filter:${tag}
