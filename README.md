**Использование**
- для компиляции и запуска: `npm start`
- только компиляция: `npm run build`
- запуск: `npm run calc`

**Установка**

Программа написана на языке TypeScript. Запускается под NodeJs.
Для запуска под Windows потребуется установка NodeJs: https://nodejs.org/en/download/

Далее необходимо установить пакеты npm.
Откройте в папке с программой терминал/консоль и выполните `npm i`

Также потребуется установить глобально TypeScript: `npm i -g typescript`

Проверьте правильность установки, узнав текущие версии установленных программ:
- `npm -v`
- `node -v`
- `tsc -v`

Для удобства работы с кодом и типами советую использовать IDE, например VSCode:
https://code.visualstudio.com/

Если возникает ошибка "Невозможно загрузить файл *.ps1, так как выполнение сценариев отключено в этой системе."
- Открываем терминал PowerShell от админа.
- Вставляем и запускаем - `Set-ExecutionPolicy RemoteSigned`
- На вопрос отвечаем - `A`
