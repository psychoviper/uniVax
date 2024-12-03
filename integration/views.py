from django.shortcuts import render
from django import forms
from django.http import HttpResponse
import os
from django.conf import settings
from .combined import f_to_c


# Create your views here.
class UploadFileForm(forms.Form):
    input_file = forms.FileField(label='Upload File', widget=forms.ClearableFileInput(attrs={'class': 'form-control'}))

def index(request):
    if request.method == "POST":
        form = UploadFileForm(request.POST)
        print('uploaded')
        # if form.is_valid():
        uploaded_file = request.FILES.get('input_file')
        print("Start")
        handle_uploaded_file(uploaded_file)
        print("End")
        return HttpResponse('File uploaded successfully!')
    else :
        form=UploadFileForm()
    return render(request, 'index.html', {'form':form})


def handle_uploaded_file(file):
    file.name='input.fa'
    print('Adding file')
    media_path = os.path.join(settings.MEDIA_ROOT)
    if not os.path.exists(media_path):
        os.makedirs(media_path)

    with open(os.path.join(media_path, file.name), 'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)
        print('Done')
    f_to_c(os.path.join(media_path, file.name))